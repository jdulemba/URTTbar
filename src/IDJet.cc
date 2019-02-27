#include "Analyses/URTTbar/interface/IDJet.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"

const std::unordered_map<std::string, IDJet::BTag> IDJet::tag_names = {
  {"NONE",      IDJet::BTag::NONE},
  {"DEEPCSVLOOSE",  IDJet::BTag::DEEPCSVLOOSE }, 
  {"DEEPCSVMEDIUM", IDJet::BTag::DEEPCSVMEDIUM},
  {"DEEPCSVTIGHT",  IDJet::BTag::DEEPCSVTIGHT }, 
  {"DEEPJETLOOSE",  IDJet::BTag::DEEPJETLOOSE }, 
  {"DEEPJETMEDIUM", IDJet::BTag::DEEPJETMEDIUM},
  {"DEEPJETTIGHT",  IDJet::BTag::DEEPJETTIGHT }, 
  {"DEEPCSVCTAGLOOSE" , IDJet::BTag::DEEPCSVCTAGLOOSE }, 
  {"DEEPCSVCTAGMEDIUM", IDJet::BTag::DEEPCSVCTAGMEDIUM}, 
  {"DEEPCSVCTAGTIGHT" , IDJet::BTag::DEEPCSVCTAGTIGHT },
  {"DEEPJETCTAGLOOSE" , IDJet::BTag::DEEPJETCTAGLOOSE }, 
  {"DEEPJETCTAGMEDIUM", IDJet::BTag::DEEPJETCTAGMEDIUM}, 
  {"DEEPJETCTAGTIGHT" , IDJet::BTag::DEEPJETCTAGTIGHT },
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
  case DEEPCSVLOOSE:
  case DEEPCSVMEDIUM:
  case DEEPCSVTIGHT: return IDType::DEEPCSV;
  case DEEPJETLOOSE:
  case DEEPJETMEDIUM:
  case DEEPJETTIGHT: return IDType::DEEPJET;
  case DEEPCSVCTAGLOOSE: 
  case DEEPCSVCTAGMEDIUM: 
  case DEEPCSVCTAGTIGHT: return IDType::DEEPCSVCTAG;
  case DEEPJETCTAGLOOSE: 
  case DEEPJETCTAGMEDIUM: 
  case DEEPJETCTAGTIGHT: return IDType::DEEPJETCTAG;
  default: return IDType::NOTSET;
  }
  return IDType::NOTSET;
}

std::string IDJet::id_string(BTag id) {
  switch(id) {
  case DEEPCSVLOOSE:
  case DEEPCSVMEDIUM:
  case DEEPCSVTIGHT: return "deepcsv";
  case DEEPJETLOOSE:
  case DEEPJETMEDIUM:
  case DEEPJETTIGHT: return "deepjet";
  case DEEPCSVCTAGLOOSE: 
  case DEEPCSVCTAGMEDIUM: 
  case DEEPCSVCTAGTIGHT: return "deepcsvctag";
  case DEEPJETCTAGLOOSE: 
  case DEEPJETCTAGMEDIUM: 
  case DEEPJETCTAGTIGHT: return "deepjetctag";
  default: return "";
  }
  return "";
}

BTagEntry::OperatingPoint IDJet::tag_tightness(BTag id) {
  BTagEntry::OperatingPoint val = BTagEntry::OperatingPoint::OP_NOTSET;
  switch(id) {
  case BTag::DEEPCSVLOOSE : val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::DEEPCSVMEDIUM: val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::DEEPCSVTIGHT : val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  case BTag::DEEPJETLOOSE : val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::DEEPJETMEDIUM: val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::DEEPJETTIGHT : val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  case BTag::DEEPCSVCTAGLOOSE:  val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::DEEPCSVCTAGMEDIUM: val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::DEEPCSVCTAGTIGHT:  val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  case BTag::DEEPJETCTAGLOOSE:  val = BTagEntry::OperatingPoint::OP_LOOSE;  break;
  case BTag::DEEPJETCTAGMEDIUM: val = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case BTag::DEEPJETCTAGTIGHT:  val = BTagEntry::OperatingPoint::OP_TIGHT;  break;
  default:  val = BTagEntry::OperatingPoint::OP_NOTSET; break;
  }
  return val;
}


bool IDJet::BTagId(BTag wp) const {
  //  cout << "working point: " << wp << endl;
  if(wp == BTag::NONE) return true;
  else if(wp == BTag::DEEPCSVLOOSE ) return DeepCSVbDisc() > 0.1241;
  else if(wp == BTag::DEEPCSVMEDIUM) return DeepCSVbDisc() > 0.4184;
  else if(wp == BTag::DEEPCSVTIGHT ) return DeepCSVbDisc() > 0.7527;
  else if(wp == BTag::DEEPJETLOOSE ) return DeepJetbDisc() > 0.0494;
  else if(wp == BTag::DEEPJETMEDIUM) return DeepJetbDisc() > 0.2770;
  else if(wp == BTag::DEEPJETTIGHT ) return DeepJetbDisc() > 0.7264;
  else {
    Logger::log().fatal() << wp << "Is not a valid b-tagging working point!"<< std::endl;
    throw 42;
  }
}


float IDJet::DeepCSVCvsLtag() const {
    return ( (DeepCSVProbC() != -1) ? (DeepCSVProbC())/(DeepCSVProbC() + DeepCSVProbUDSG()) : -1 );
}

float IDJet::DeepCSVCvsBtag() const {
    return ( (DeepCSVProbC() != -1) ? (DeepCSVProbC())/(DeepCSVProbC() + DeepCSVProbB() + DeepCSVProbBB()) : -1 );
}

float IDJet::DeepJetCvsLtag() const {
    return ( (DeepFlavourProbC() != -1) ? (DeepFlavourProbC())/(DeepFlavourProbC() + DeepFlavourProbUDS() + DeepFlavourProbG()) : -1 );
}

float IDJet::DeepJetCvsBtag() const {
    return ( (DeepFlavourProbC() != -1) ? (DeepFlavourProbC())/(DeepFlavourProbC() + DeepFlavourProbB() + DeepFlavourProbBB() + DeepFlavourProbLepB()) : -1 );
}

bool IDJet::CTagId(BTag wp) const	{
  double cvsl_thr = -1.;
  double cvsb_thr = -1.;
  if(wp == BTag::NONE) return true;
  else if(wp == BTag::DEEPCSVCTAGLOOSE)  {cvsl_thr = 0.040; cvsb_thr = 0.35; return (DeepCSVCvsLtag() > cvsl_thr && DeepCSVCvsBtag() > cvsb_thr);}
  else if(wp == BTag::DEEPCSVCTAGMEDIUM) {cvsl_thr = 0.137; cvsb_thr = 0.29; return (DeepCSVCvsLtag() > cvsl_thr && DeepCSVCvsBtag() > cvsb_thr);}
  else if(wp == BTag::DEEPCSVCTAGTIGHT)  {cvsl_thr = 0.660; cvsb_thr = 0.10; return (DeepCSVCvsLtag() > cvsl_thr && DeepCSVCvsBtag() > cvsb_thr);}
  else if(wp == BTag::DEEPJETCTAGLOOSE)  {cvsl_thr = 0.030; cvsb_thr = 0.40; return (DeepJetCvsLtag() > cvsl_thr && DeepJetCvsBtag() > cvsb_thr);}
  else if(wp == BTag::DEEPJETCTAGMEDIUM) {cvsl_thr = 0.085; cvsb_thr = 0.29; return (DeepJetCvsLtag() > cvsl_thr && DeepJetCvsBtag() > cvsb_thr);}
  else if(wp == BTag::DEEPJETCTAGTIGHT)  {cvsl_thr = 0.480; cvsb_thr = 0.05; return (DeepJetCvsLtag() > cvsl_thr && DeepJetCvsBtag() > cvsb_thr);}
  else {
    Logger::log().fatal() << wp << "Is not a valid C-tagging working point!"<< std::endl;
    throw 42;
  }
}

bool IDJet::TagId(BTag wp) const {
  switch(id_type(wp)) {
  case DEEPCSV: return BTagId(wp);
  case DEEPJET: return BTagId(wp);
  case DEEPCSVCTAG: return CTagId(wp);
  case DEEPJETCTAG: return CTagId(wp);
  default: return false;
  }
}
