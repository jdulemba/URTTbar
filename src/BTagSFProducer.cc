#include "Analyses/URTTbar/interface/BTagSFProducer.h"
#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "URAnalysis/AnalysisFW/interface/PDGID.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include <string>
#include "TFile.h"
#include <cmath>
using namespace systematics;

BTagSFProducer::BTagSFProducer(TTPermutator &permutator, float float_c, float float_l, float float_b):
  eff_light_loose(), 
  eff_light_tight(),
  eff_charm_loose(), 
  eff_charm_tight(),
  eff_bottom_loose(), 
  eff_bottom_tight()
{
  URParser &parser = URParser::instance();
  parser.addCfgParameter<std::string>("general", "btag_sf", "source file for btag scale factors");
  parser.addCfgParameter<std::string>("general", "btag_eff", "source file for btag efficiencies");
  parser.parseArguments();

  DataFile sf_file( parser.getCfgPar<std::string>("general", "btag_sf"));
  DataFile eff_file(parser.getCfgPar<std::string>("general", "btag_eff"));
  configure(sf_file, eff_file, permutator.tight_bID_cut(), permutator.loose_bID_cut(), float_c, float_l, float_b);
}

BTagSFProducer::BTagSFProducer(std::string tight, std::string loose, float float_c, float float_l, float float_b): 
  eff_light_loose(), 
  eff_light_tight(),
  eff_charm_loose(), 
  eff_charm_tight(),
  eff_bottom_loose(), 
  eff_bottom_tight()
{
  URParser &parser = URParser::instance();
  parser.addCfgParameter<std::string>("general", "btag_sf", "source file for btag scale factors");
  parser.addCfgParameter<std::string>("general", "btag_eff", "source file for btag efficiencies");

  size_t dot = tight.find(".");
  std::string group(tight, 0, dot);
  std::string par(tight, dot+1, tight.size()-dot);
  parser.addCfgParameter<std::string>(group, par, "");
  if(loose.size() > 0) {
    dot = loose.find(".");
    std::string group2(loose, 0, dot);
    std::string par2(loose, dot+1, loose.size()-dot);
    parser.addCfgParameter<std::string>(group2, par2, "");
  }

  parser.parseArguments();

  string a = parser.getCfgPar<string>(tight);
  IDJet::BTag tighttag = IDJet::tag(a);

  IDJet::BTag loosetag = IDJet::BTag::NONE;
  if(loose.size() > 0) {
    string b = parser.getCfgPar<string>(loose);
    //cout << loose << " --> " << b << endl;
    loosetag = IDJet::tag(b);
  }

  DataFile sf_file( parser.getCfgPar<std::string>("general", "btag_sf"));
  DataFile eff_file(parser.getCfgPar<std::string>("general", "btag_eff"));
  configure(sf_file, eff_file, tighttag, loosetag, float_c, float_l, float_b);
}

BTagSFProducer::BTagSFProducer(const DataFile &sf_file, const DataFile &eff_file,
                               IDJet::BTag tighttag, IDJet::BTag loosetag, 
                               float float_c, float float_l, float float_b): 
  eff_light_loose(), 
  eff_light_tight(),
  eff_charm_loose(), 
  eff_charm_tight(),
  eff_bottom_loose(), 
  eff_bottom_tight()
{
	Logger::log().debug() << "BTagSFProducer ctor" << endl;
  configure(sf_file, eff_file, tighttag, loosetag, float_c, float_l, float_b);
}

void BTagSFProducer::configure(const DataFile &sf_file, const DataFile &eff_file,
                               IDJet::BTag tighttag, IDJet::BTag loosetag, 
                               float float_c, float float_l, float float_b) {
  float_c_ = float_c;
  float_l_ = float_l;
  float_b_ = float_b;
  tight_ = tighttag;
  loose_ = loosetag;

  if(float_c_ > 1. || float_l_ > 1. || float_b_ > 1.) {
    Logger::log().fatal() << "You cannot float a SF for a value > 1!" << std::endl;
    throw 42;
  }

  //Get WP used
  //defined by the permutator they MUST be the same, therefore, not double definition
  BTagEntry::OperatingPoint wp_tight = IDJet::tag_tightness(tighttag);  
  BTagEntry::OperatingPoint wp_loose = IDJet::tag_tightness(loosetag);
  no_loose_cut_ = (wp_loose == BTagEntry::OperatingPoint::OP_NOTSET);
  if(no_loose_cut_) wp_loose = wp_tight;

  //Get efficiencies
  TFile eff_tfile(eff_file.path().c_str());
  TH1::AddDirectory(false);
  eff_light_tight  = get_from<TH2D>(eff_tfile, "light/"  + IDJet::tag2string(tighttag) + "_eff", "t_l_eff");
  eff_charm_tight  = get_from<TH2D>(eff_tfile, "charm/"  + IDJet::tag2string(tighttag) + "_eff", "t_c_eff");
  eff_bottom_tight = get_from<TH2D>(eff_tfile, "bottom/" + IDJet::tag2string(tighttag) + "_eff", "t_b_eff");
  
  if(!no_loose_cut_){
    eff_light_loose  = get_from<TH2D>(eff_tfile, "light/"  + IDJet::tag2string(loosetag) + "_eff", "l_l_eff");
    eff_charm_loose  = get_from<TH2D>(eff_tfile, "charm/"  + IDJet::tag2string(loosetag) + "_eff", "l_c_eff");
    eff_bottom_loose = get_from<TH2D>(eff_tfile, "bottom/" + IDJet::tag2string(loosetag) + "_eff", "l_b_eff");
  }
  else {
    eff_light_loose  = eff_light_tight ;
    eff_charm_loose  = eff_charm_tight ;
    eff_bottom_loose = eff_bottom_tight;
  }    
  eff_tfile.Close();
  TH1::AddDirectory(true);

  if(float_c_ >= 0 && float_l_ >= 0 && float_b_ >= 0) return; //No SF to be used, so just skip, they might not even be present!

  //Get SFs
  std::string calibration_type = IDJet::id_string(tighttag);
  if(!no_loose_cut_) {
    if(calibration_type != IDJet::id_string(loosetag)) {
      Logger::log().fatal() << "BTagSFProducer: Mixing cuts on " << calibration_type << 
        "and " << IDJet::id_string(loosetag) << " is not allowed!" << std::endl;
      throw 42;
    }
  }

  //
  try {
    calibration_ = BTagCalibration(calibration_type, sf_file.path());
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
}

BTagSFProducer::~BTagSFProducer() {
}

double BTagSFProducer::scale_factor(const std::vector<IDJet*> &jets, systematics::SysShifts shift) {
  if(ignore_partial_shifts_ && 
     (shift == SysShifts::BTAG_B_UP || shift == SysShifts::BTAG_B_DW || shift == SysShifts::BTAG_C_UP || 
      shift == SysShifts::BTAG_C_DW || shift == SysShifts::BTAG_L_UP || shift == SysShifts::BTAG_L_DW)) {
    shift = SysShifts::NOSYS; //reset value
  }
  if(ignore_general_shifts_ && (shift == SysShifts::BTAG_UP || shift == SysShifts::BTAG_DW))  {
    shift = SysShifts::NOSYS;
  }

  double mc_prob=1;
  double data_like_prob=1; //it's called data, but is on MC!

  BTagCalibrationReader *reader_loose = 0; //&readers_loose_[1];
  BTagCalibrationReader *reader_tight = 0; //&readers_tight_[1];

  for(auto jet : jets) {
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
      reader_loose = &readers_loose_[idx];
      reader_tight = &readers_tight_[idx];
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

    if(jet->TagId(tight_)) {
      mc_prob *= eff_tight;
      data_like_prob *= eff_tight*tight_sf;
    }
    else if(loose_ != IDJet::BTag::NONE && jet->TagId(loose_)) {
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
