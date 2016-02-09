#ifndef systematics_h
#define systematics_h

#include <map>
#include <string>
#include <vector>

namespace systematics {
  enum SysShifts {
    NOSYS, JES_UP, JES_DW, JER_UP, JER_DW, MET_UP, MET_DW, 
    RENORM_UP, RENORM_DW, FACTOR_UP, FACTOR_DW,
    BTAG_UP, BTAG_DW, BTAG_B_UP, BTAG_B_DW, BTAG_C_UP, BTAG_C_DW, BTAG_L_UP, BTAG_L_DW
  };
  const std::map<std::string, SysShifts> name_to_shift = {
    {"nosys", SysShifts::NOSYS}, 
    {"jes_up", SysShifts::JES_UP}, 
    {"jes_down", SysShifts::JES_DW}, 
    {"jer_up", SysShifts::JER_UP},
    {"jer_down", SysShifts::JER_DW}, 
    {"met_up", SysShifts::MET_UP},
    {"met_down", SysShifts::MET_DW},
    {"renorm_up", SysShifts::RENORM_UP}, 
    {"renorm_down", SysShifts::RENORM_DW}, 
    {"factor_up", SysShifts::FACTOR_UP}, 
    {"factor_down", SysShifts::FACTOR_DW},
    {"btag_up", SysShifts::BTAG_UP}, 
    {"btag_down", SysShifts::BTAG_DW}, 
    {"btagb_up", SysShifts::BTAG_B_UP}, 
    {"btagb_down", SysShifts::BTAG_B_DW}, 
    {"btagc_up", SysShifts::BTAG_C_UP}, 
    {"btagc_down", SysShifts::BTAG_C_DW}, 
    {"btagl_up", SysShifts::BTAG_L_UP}, 
    {"btagl_down", SysShifts::BTAG_L_DW}
  };
  const std::map<SysShifts, std::string> shift_to_name = {
    {SysShifts::NOSYS , "nosys"   }, 
    {SysShifts::JES_UP, "jes_up"  }, 
    {SysShifts::JES_DW, "jes_down"}, 
    {SysShifts::JER_UP, "jer_up"  },
    {SysShifts::JER_DW, "jer_down"}, 
    {SysShifts::MET_UP, "met_up"  },
    {SysShifts::MET_DW, "met_down"},
    {SysShifts::RENORM_UP, "renorm_up"  }, 
    {SysShifts::RENORM_DW, "renorm_down"}, 
    {SysShifts::FACTOR_UP, "factor_up"  }, 
    {SysShifts::FACTOR_DW, "factor_down"},
    {SysShifts::BTAG_UP  , "btag_up"   }, 
    {SysShifts::BTAG_DW  , "btag_down" }, 
    {SysShifts::BTAG_B_UP, "btagb_up"  }, 
    {SysShifts::BTAG_B_DW, "btagb_down"}, 
    {SysShifts::BTAG_C_UP, "btagc_up"  }, 
    {SysShifts::BTAG_C_DW, "btagc_down"}, 
    {SysShifts::BTAG_L_UP, "btagl_up"  }, 
    {SysShifts::BTAG_L_DW, "btagl_down"}
  };

  std::vector<SysShifts> get_systematics(std::string outname);
  std::string get_sample(std::string outname);
};

#endif
