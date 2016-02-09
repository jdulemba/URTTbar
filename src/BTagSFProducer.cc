#include "BTagSFProducer.h"
#include "TTPermutator.h"
#include "PDGID.h"
#include "URParser.h"
#include <string>
#include "DataFile.h"
#include "Logger.h"
#include "TFile.h"
#include <cmath>

BTagSFProducer::BTagSFProducer(TTPermutator &permutator, float float_c, float float_l, float float_b):
  float_c_(float_c),
  float_l_(float_l),
  float_b_(float_b) {

  if(float_c_ > 1. || float_l_ > 1. || float_b_ > 1.) {
    Logger::log().fatal() << "You cannot float a SF for a value > 1!" << std::endl;
    throw 42;
  }

  URParser &parser = URParser::instance();
  parser.addCfgParameter<std::string>("general", "btag_sf", "source file for btag scale factors");
  parser.addCfgParameter<std::string>("general", "btag_eff", "source file for btag efficiencies");
  parser.parseArguments();

  DataFile sf_file(parser.getCfgPar<std::string>("general", "btag_sf"));
  DataFile eff_file(parser.getCfgPar<std::string>("general", "btag_eff"));

  //defined by the permutator they MUST be the same, therefore, not double definition
  BTagEntry::OperatingPoint wp_tight;  
  switch(permutator.tight_bID_cut()){
  case IDJet::BTag::CSVLOOSE: wp_tight = BTagEntry::OperatingPoint::OP_LOOSE; break;
  case IDJet::BTag::CSVMEDIUM: wp_tight = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case IDJet::BTag::CSVTIGHT: wp_tight = BTagEntry::OperatingPoint::OP_TIGHT; break;
  case IDJet::BTag::NONE: break;
  }

  BTagEntry::OperatingPoint wp_loose;
  switch(permutator.loose_bID_cut()){
  case IDJet::BTag::CSVLOOSE: wp_loose = BTagEntry::OperatingPoint::OP_LOOSE; break;
  case IDJet::BTag::CSVMEDIUM: wp_loose = BTagEntry::OperatingPoint::OP_MEDIUM; break;
  case IDJet::BTag::CSVTIGHT: wp_loose = BTagEntry::OperatingPoint::OP_TIGHT; break;
  case IDJet::BTag::NONE: wp_loose = wp_tight; no_loose_cut_=true; break;
  }

  //
  try {
    calibration_ = BTagCalibration("csvv2", sf_file.path());
  } catch(std::exception e) {
    Logger::log().fatal() << "BTagSFProducer::BTagSFProducer caught an exception in instantiating the BTagCalibration with input file: " <<
      sf_file.path() << std::endl;
    throw 42;
  }
  readers_tight_[0] = BTagCalibrationReader(&calibration_, wp_tight, "mujets", "down"); //[down, central, up]
  readers_tight_[1] = BTagCalibrationReader(&calibration_, wp_tight, "mujets", "central"); //[down, central, up]
  readers_tight_[2] = BTagCalibrationReader(&calibration_, wp_tight, "mujets", "up"); //[down, central, up]

  readers_loose_[0] = BTagCalibrationReader(&calibration_, wp_loose, "mujets", "down"); //[down, central, up]
  readers_loose_[1] = BTagCalibrationReader(&calibration_, wp_loose, "mujets", "central"); //[down, central, up]
  readers_loose_[2] = BTagCalibrationReader(&calibration_, wp_loose, "mujets", "up"); //[down, central, up]

  TH1::AddDirectory(false);
  TFile eff_tfile(eff_file.path().c_str());
  eff_light_tight  = get_from<TH2D>(eff_tfile, "light/tight_eff", "t_eff");
  eff_charm_tight  = get_from<TH2D>(eff_tfile, "charm/tight_eff", "t_eff");
  eff_bottom_tight = get_from<TH2D>(eff_tfile, "bottom/tight_eff", "t_eff");
  
  if(!no_loose_cut_){
    eff_light_loose = get_from<TH2D>(eff_tfile, "light/loose_eff", "l_eff");
    eff_charm_loose = get_from<TH2D>(eff_tfile, "charm/loose_eff", "l_eff");
    eff_bottom_loose = get_from<TH2D>(eff_tfile, "bottom/loose_eff", "l_eff");
  }
  else {
    eff_light_loose  = eff_light_tight ;
    eff_charm_loose  = eff_charm_tight ;
    eff_bottom_loose = eff_bottom_tight;
  }    
  TH1::AddDirectory(true);
}

BTagSFProducer::~BTagSFProducer() {
  if(eff_light_tight) {
    delete eff_light_tight ;
    delete eff_charm_tight ;
    delete eff_bottom_tight;
    if(!no_loose_cut_) {
      delete eff_light_loose ;
      delete eff_charm_loose ;
      delete eff_bottom_loose;  
    }    
  }
}

double BTagSFProducer::scale_factor(TTPermutator &permutator, systematics::SysShifts shift) {
  double mc_prob=1;
  double data_like_prob=1; //it's called data, but is on MC!

  BTagCalibrationReader *reader_loose = 0; //&readers_loose_[1];
  BTagCalibrationReader *reader_tight = 0; //&readers_tight_[1];

  auto loose_b = permutator.loose_bID_cut();
  auto tight_b = permutator.tight_bID_cut();
  
  for(auto jet : permutator.capped_jets()) {
    BTagEntry::JetFlavor jet_flav;
    double jpt = jet->Pt();
    double jeta = jet->Eta();
    double eff_loose=0;
    double eff_tight=0;

    //ASSUME ALL EFF HISTOS HAVE SAME BINNING!
    int binx = eff_bottom_loose->GetXaxis()->FindFixBin(jpt);
    binx = TMath::Min(binx, eff_bottom_loose->GetNbinsX());
    int biny = eff_bottom_loose->GetYaxis()->FindFixBin(jeta);
    biny = TMath::Min(biny, eff_bottom_loose->GetNbinsY());
    float float_value=0;

    switch(TMath::Abs(jet->hadronFlavour())) {
    case ura::PDGID::b: 
      if(float_b_ == 0) continue;
      float_value = float_b_;
      jet_flav=BTagEntry::JetFlavor::FLAV_B;
      eff_loose = eff_bottom_loose->GetBinContent(binx, biny);
      eff_tight = eff_bottom_tight->GetBinContent(binx, biny);
      break;
    case ura::PDGID::c: 
      if(float_c_ == 0) continue;
      float_value = float_c_;
      jet_flav=BTagEntry::JetFlavor::FLAV_C;       
      eff_loose = eff_charm_loose->GetBinContent(binx, biny);
      eff_tight = eff_charm_tight->GetBinContent(binx, biny);
      break;
    default: 
      if(float_l_ == 0) continue;
      float_value = float_l_;
      jet_flav=BTagEntry::JetFlavor::FLAV_UDSG; 
      eff_loose = eff_light_loose->GetBinContent(binx, biny);
      eff_tight = eff_light_tight->GetBinContent(binx, biny);
      break;
    } //switch(Abs(jet->hadronFlavour()))  

    unsigned int idx=1;
    if(shift == systematics::SysShifts::BTAG_UP) idx=2;
    else if(shift == systematics::SysShifts::BTAG_DW) idx=0;
    else if(shift == systematics::SysShifts::BTAG_L_UP && jet_flav == BTagEntry::JetFlavor::FLAV_UDSG) idx=2;
    else if(shift == systematics::SysShifts::BTAG_B_UP && jet_flav == BTagEntry::JetFlavor::FLAV_B) idx=2;
    else if(shift == systematics::SysShifts::BTAG_C_UP && jet_flav == BTagEntry::JetFlavor::FLAV_C) idx=2;
    else if(shift == systematics::SysShifts::BTAG_L_DW && jet_flav == BTagEntry::JetFlavor::FLAV_UDSG) idx=0;
    else if(shift == systematics::SysShifts::BTAG_B_DW && jet_flav == BTagEntry::JetFlavor::FLAV_B) idx=0;
    else if(shift == systematics::SysShifts::BTAG_C_DW && jet_flav == BTagEntry::JetFlavor::FLAV_C) idx=0;

    // if( 
    //   (skip_b_ && jet_flav == BTagEntry::JetFlavor::FLAV_B) ||
    //   (skip_c_ && jet_flav == BTagEntry::JetFlavor::FLAV_C) || 
    //   (skip_l_ && jet_flav == BTagEntry::JetFlavor::FLAV_UDSG)
    //   ) continue;

    //access SF now, so we can catch exceptions!
    double tight_sf = -1.;
    double loose_sf = -1.;

    if(float_value < 0) { //use provided SF
      BTagCalibrationReader *reader_loose = &readers_loose_[idx];
      BTagCalibrationReader *reader_tight = &readers_tight_[idx];
      try { 
        tight_sf = reader_tight->eval(jet_flav, jet->Eta(), jet->Pt());
        loose_sf = reader_loose->eval(jet_flav, jet->Eta(), jet->Pt());
      } catch(std::out_of_range e) {
        Logger::log().fatal() << "Problem accessing BTV SF for jet: " << jet_flav <<
          ", " << jet->Eta() << ", " << jet->Pt() << std::endl;
        throw 42;
      }
    }
    else { //use forced value
      if(idx == 0) { //down
        tight_sf = 1-float_value;
        loose_sf = tight_sf;
      }
      else if(idx == 1) {//center
        tight_sf = 1.;
        loose_sf = 1.;
      }
      else { //up
        tight_sf = 1+float_value;
        loose_sf = 1+float_value;
      }
    }

    if(jet->BTagId(tight_b)) {
      mc_prob *= eff_tight;
      data_like_prob *= eff_tight*tight_sf;
    }
    else if(loose_b != IDJet::BTag::NONE && jet->BTagId(loose_b)) {
      mc_prob *= (eff_loose - eff_tight);
      data_like_prob *= (
        eff_loose*loose_sf - 
        eff_tight*tight_sf
        );
    }
    else {
      mc_prob *= (1-eff_loose);
      data_like_prob *= (1 - eff_loose*loose_sf);
    }
  } //for(auto jet : permutator.capped_jets())
  
  if(mc_prob == 0) {
    Logger::log().error() << "MC Probability is 0!" << std::endl;
    throw 49;
  }
  return (mc_prob != 0.) ? data_like_prob/mc_prob : 0.;
}
