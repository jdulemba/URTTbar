#include "Analyses/URTTbar/interface/IDJet.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"

const std::unordered_map<std::string, IDJet::BTag> IDJet::tag_names = {
  {"NONE",      IDJet::BTag::NONE},
  {"CSVLOOSE",  IDJet::BTag::CSVLOOSE},
  {"CSVMEDIUM", IDJet::BTag::CSVMEDIUM},
  {"CSVTIGHT",  IDJet::BTag::CSVTIGHT},
  {"MVALOOSE",  IDJet::BTag::MVALOOSE},
  {"MVAMEDIUM", IDJet::BTag::MVAMEDIUM},
  {"MVATIGHT",  IDJet::BTag::MVATIGHT},
  {"CTAGLOOSE" , IDJet::BTag::CTAGLOOSE }, 
  {"CTAGMEDIUM", IDJet::BTag::CTAGMEDIUM}, 
  {"CTAGTIGHT" , IDJet::BTag::CTAGTIGHT },
  //{"DEEPCTAGLOOSE" , IDJet::BTag::DEEPCTAGLOOSE }, 
  //{"DEEPCTAGMEDIUM", IDJet::BTag::DEEPCTAGMEDIUM}, 
  //{"DEEPCTAGTIGHT" , IDJet::BTag::DEEPCTAGTIGHT },
  {"DEEPCSVLOOSE",  IDJet::BTag::DEEPCSVLOOSE }, 
  {"DEEPCSVMEDIUM", IDJet::BTag::DEEPCSVMEDIUM},
  {"DEEPCSVTIGHT",  IDJet::BTag::DEEPCSVTIGHT }, 
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
  case MVALOOSE: 
  case MVAMEDIUM: 
  case MVATIGHT: return IDType::MVA;
  case CTAGLOOSE: 
  case CTAGMEDIUM: 
  case CTAGTIGHT: return IDType::CTAG;
  //case DEEPCTAGLOOSE: 
  //case DEEPCTAGMEDIUM: 
  //case DEEPCTAGTIGHT: return IDType::DEEPCTAG;
  case DEEPCSVLOOSE:
  case DEEPCSVMEDIUM:
  case DEEPCSVTIGHT: return IDType::DEEPCSV;
  default: return IDType::NOTSET;
  }
  return IDType::NOTSET;
}

std::string IDJet::id_string(BTag id) {
  switch(id) {
  case CSVLOOSE: 
  case CSVMEDIUM: 
  case CSVTIGHT: return "csvv2";
  case MVALOOSE: 
  case MVAMEDIUM: 
  case MVATIGHT: return "cMVAv2";
  case CTAGLOOSE: 
  case CTAGMEDIUM: 
  case CTAGTIGHT: return "ctag";
  //case DEEPCTAGLOOSE: 
  //case DEEPCTAGMEDIUM: 
  //case DEEPCTAGTIGHT: return "deepctag";
  case DEEPCSVLOOSE:
  case DEEPCSVMEDIUM:
  case DEEPCSVTIGHT: return "deepcsv";
  default: return "";
  }
  return "";
}

BTagEntry::OperatingPoint IDJet::tag_tightness(BTag id) {
  BTagEntry::OperatingPoint val = BTagEntry::OperatingPoint::OP_NOTSET;
  switch(id) {
  case BTag::CSVLOOSE:   val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::CTAGLOOSE:  val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::MVALOOSE:   val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::MVAMEDIUM:  val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::CTAGMEDIUM: val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::CSVMEDIUM:  val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::CSVTIGHT:   val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  case BTag::CTAGTIGHT:  val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  case BTag::MVATIGHT:   val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  //case BTag::DEEPCTAGLOOSE:  val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  //case BTag::DEEPCTAGMEDIUM: val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  //case BTag::DEEPCTAGTIGHT:  val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  case BTag::DEEPCSVLOOSE : val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::DEEPCSVMEDIUM: val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::DEEPCSVTIGHT : val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  default:  val = BTagEntry::OperatingPoint::OP_NOTSET; break;
  }
  return val;
}


bool IDJet::BTagId(BTag wp) const {
  //  cout << "working point: " << wp << endl;
  if(wp == BTag::NONE) return true;
  //else if(wp == BTag::CSVLOOSE)  return csvIncl() > 0.5803;
  //else if(wp == BTag::CSVMEDIUM) return csvIncl() > 0.8838;
  //else if(wp == BTag::CSVTIGHT)  return csvIncl() > 0.9693;
  //else if(wp == BTag::MVALOOSE)  return CombinedMVA() > -0.5884;
  //else if(wp == BTag::MVAMEDIUM) return CombinedMVA() > 0.4432;
  //else if(wp == BTag::MVATIGHT)  return CombinedMVA() > 0.9432;
  //else if(wp == BTag::DEEPCSVLOOSE ) return (DeepCSVProbB() + DeepCSVProbBB()) > 0.1522;
  //else if(wp == BTag::DEEPCSVMEDIUM) return (DeepCSVProbB() + DeepCSVProbBB()) > 0.4941;
  //else if(wp == BTag::DEEPCSVTIGHT ) return (DeepCSVProbB() + DeepCSVProbBB()) > 0.8001;
  else if(wp == BTag::CSVLOOSE)  return csvIncl() > 0.5426;
  else if(wp == BTag::CSVMEDIUM) return csvIncl() > 0.8484;
  else if(wp == BTag::CSVTIGHT)  return csvIncl() > 0.9535;
  else if(wp == BTag::MVALOOSE)  return CombinedMVA() > -0.5884;
  else if(wp == BTag::MVAMEDIUM) return CombinedMVA() > 0.4432;
  else if(wp == BTag::MVATIGHT)  return CombinedMVA() > 0.9432;
  else if(wp == BTag::DEEPCSVLOOSE ) return (DeepCSVProbB() + DeepCSVProbBB()) > 0.2219;
  else if(wp == BTag::DEEPCSVMEDIUM) return (DeepCSVProbB() + DeepCSVProbBB()) > 0.6324;
  else if(wp == BTag::DEEPCSVTIGHT ) return (DeepCSVProbB() + DeepCSVProbBB()) > 0.8958;
  else {
    Logger::log().fatal() << wp << "Is not a valid b-tagging working point!"<< std::endl;
    throw 42;
  }
}


//float IDJet::DeepCSVCvsLtag() const {
//    return ( (DeepCSVProbC() != -1) ? (DeepCSVProbC())/(DeepCSVProbC() + DeepCSVProbUDSG()) : -1 );
//}
//
//float IDJet::DeepCSVCvsBtag() const {
//    return ( (DeepCSVProbC() != -1) ? (DeepCSVProbC())/(DeepCSVProbC() + DeepCSVProbB() + DeepCSVProbBB()) : -1 );
//}

bool IDJet::CTagId(BTag wp) const	{
  double cvsl_thr = -1.;
  double cvsb_thr = -1.;
  //cout << "Deep: " << DeepCSVCvsBtag() << endl;
  //cout << "DeepCSV CvsL disc: " << (DeepCSVProbC())/(DeepCSVProbC() + DeepCSVProbB() + DeepCSVProbBB()) << endl;
  if(wp == BTag::NONE) return true;
  else if(wp == BTag::CTAGLOOSE)  {cvsl_thr = -0.337; cvsb_thr = -0.356; return (CvsLtag() > cvsl_thr && CvsBtag() > cvsb_thr);}
  else if(wp == BTag::CTAGMEDIUM) {cvsl_thr = -0.073; cvsb_thr = -0.302; return (CvsLtag() > cvsl_thr && CvsBtag() > cvsb_thr);}
  else if(wp == BTag::CTAGTIGHT)  {cvsl_thr = 0.294; cvsb_thr = -0.682; return (CvsLtag() > cvsl_thr && CvsBtag() > cvsb_thr);}
  //else if(wp == BTag::CTAGLOOSE)  {cvsl_thr = -0.53; cvsb_thr = -0.26; return (CvsLtag() > cvsl_thr && CvsBtag() > cvsb_thr);}
  //else if(wp == BTag::CTAGMEDIUM) {cvsl_thr = 0.07; cvsb_thr = -0.10; return (CvsLtag() > cvsl_thr && CvsBtag() > cvsb_thr);}
  //else if(wp == BTag::CTAGTIGHT)  {cvsl_thr = 0.87; cvsb_thr = -0.30; return (CvsLtag() > cvsl_thr && CvsBtag() > cvsb_thr);}
  //else if(wp == BTag::DEEPCTAGLOOSE)  {cvsl_thr = 0.05; cvsb_thr = 0.33; return (DeepCSVCvsLtag() > cvsl_thr && DeepCSVCvsBtag() > cvsb_thr);}
  //else if(wp == BTag::DEEPCTAGMEDIUM) {cvsl_thr = 0.15; cvsb_thr = 0.28; return (DeepCSVCvsLtag() > cvsl_thr && DeepCSVCvsBtag() > cvsb_thr);}
  //else if(wp == BTag::DEEPCTAGTIGHT)  {cvsl_thr = 0.80; cvsb_thr = 0.10; return (DeepCSVCvsLtag() > cvsl_thr && DeepCSVCvsBtag() > cvsb_thr);}
  else {
    Logger::log().fatal() << wp << "Is not a valid C-tagging working point!"<< std::endl;
    throw 42;
  }
  //return (CvsLtag() > cvsl_thr && CvsBtag() > cvsb_thr);
}

bool IDJet::TagId(BTag wp) const {
  switch(id_type(wp)) {
  case MVA: return BTagId(wp);
  case CSV: return BTagId(wp);
  case DEEPCSV: return BTagId(wp);
  case CTAG: return CTagId(wp);
  //case DEEPCTAG: return CTagId(wp);
  default: return false;
  }
}
