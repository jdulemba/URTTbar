#include "BTagSFProducer.h"
#include "TTPermutator.h"
#include "PDGID.h"
#include "URParser.h"
#include <string>
#include "DataFile.h"
#include "Logger.h"
#include "TFile.h"
#include <cmath>

BTagSFProducer::BTagSFProducer(TTPermutator &permutator) {
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

  //TODO to implement sys shifts
  BTagCalibrationReader *reader_loose = &readers_loose_[1];
  BTagCalibrationReader *reader_tight = &readers_tight_[1];

  auto loose_b = permutator.loose_bID_cut();
  auto tight_b = permutator.tight_bID_cut();
  
  for(auto jet : permutator.capped_jets()) {
    BTagEntry::JetFlavor jet_flav;
    double jpt = jet->Pt();
    double jeta = jet->Eta();
    double eff_loose=0;
    double eff_tight=0;
    switch(TMath::Abs(jet->hadronFlavour())) {
    case ura::PDGID::b: 
      jet_flav=BTagEntry::JetFlavor::FLAV_B;
      eff_loose = eff_bottom_loose->GetBinContent(eff_bottom_loose->FindFixBin(jpt, jeta));
      eff_tight = eff_bottom_tight->GetBinContent(eff_bottom_tight->FindFixBin(jpt, jeta));
      break;
    case ura::PDGID::c: 
      jet_flav=BTagEntry::JetFlavor::FLAV_C; 
      eff_loose = eff_charm_loose->GetBinContent(eff_charm_loose->FindFixBin(jpt, jeta));
      eff_tight = eff_charm_tight->GetBinContent(eff_charm_tight->FindFixBin(jpt, jeta));
      break;
    default: 
      jet_flav=BTagEntry::JetFlavor::FLAV_UDSG; 
      eff_loose = eff_light_loose->GetBinContent(eff_light_loose->FindFixBin(jpt, jeta));
      eff_tight = eff_light_tight->GetBinContent(eff_light_tight->FindFixBin(jpt, jeta));
      break;
    } //switch(Abs(jet->hadronFlavour()))  

    //access SF now, so we can catch exceptions!
    double tight_sf = -1.;
    double loose_sf = -1.;
    try { 
      tight_sf = reader_tight->eval(jet_flav, jet->Eta(), jet->Pt());
      loose_sf = reader_loose->eval(jet_flav, jet->Eta(), jet->Pt());
    } catch(std::out_of_range e) {
      Logger::log().fatal() << "Problem accessing BTV SF for jet: " << jet_flav <<
        ", " << jet->Eta() << ", " << jet->Pt() << std::endl;
      throw 42;
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
