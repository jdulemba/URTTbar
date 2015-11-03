#include "IDJet.h"
#include "Logger.h"

const std::unordered_map<std::string, IDJet::BTag> IDJet::tag_names = {
  {"NONE",      IDJet::BTag::NONE},
  {"CSVLOOSE",  IDJet::BTag::CSVLOOSE},
  {"CSVMEDIUM", IDJet::BTag::CSVMEDIUM},
  {"CSVTIGHT",  IDJet::BTag::CSVTIGHT},
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
