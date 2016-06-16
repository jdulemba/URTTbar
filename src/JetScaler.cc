#include "Analyses/URTTbar/interface/JetScaler.h"
#include <string>
#include <iostream>
#include <sstream>
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "URAnalysis/AnalysisFW/interface/Wards.h"

using namespace std;

JetScaler::JetScaler():
  Heta(0),
  HptsP(),
  HptsM(),
  sfc_(),
  sfu_(),
  sfd_() {
  //get configuration
  URParser &parser = URParser::instance();
  parser.addCfgParameter<string>("JERC", "JES", "");
  parser.addCfgParameter<string>("JERC", "JER", "");
    
  parser.parseArguments();

  DataFile jerf(parser.getCfgPar<string>("JERC", "JER"));
  DataFile jesf(parser.getCfgPar<string>("JERC", "JES"));

	HistoOwnershipWard hward;
	DirectoryWard dward;
  //Get JES
  TFile fjes(jesf.path().c_str());
  Heta = get_from<TH1D>(fjes, "eta", "heta"); //dynamic_cast<TH1D*>(fjes_->Get("eta"));
  for(int i = 0 ; i < Heta->GetNbinsX() ; ++i) {
    stringstream hn;
    hn << "down_" << i;
		stringstream hnn;
		hnn << "down_" << i << "_clone";
    HptsP.push_back(
			get_from<TH1D>(fjes, hn.str(), hnn.str())
			);
  }	

  for(int i = 0 ; i < Heta->GetNbinsX() ; ++i) {
    stringstream hn;
    hn << "up_" << i;
		stringstream hnn;
		hnn << "up_" << i << "_clone";
    HptsM.push_back(
			get_from<TH1D>(fjes, hn.str(), hnn.str())
			);
  }	

  //Get JER
  TFile fjer(jerf.path().c_str());
  res_ = get_from<TF1 >(fjer, "resolution", "clone_resolution");
  sfc_ = get_from<TH1F>(fjer, "sf_central", "clone_sf_central");
  sfu_ = get_from<TH1F>(fjer, "sf_up"     , "clone_sf_up"     );
  sfd_ = get_from<TH1F>(fjer, "sf_down"   , "clone_sf_down"   );
}
