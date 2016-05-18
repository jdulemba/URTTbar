#include "Analyses/URTTbar/interface/JetScaler.h"
#include <string>
#include <iostream>
#include <sstream>
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"

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

  //Get JES
  TDirectory* olddir = gDirectory;
  fjes_ = new TFile(jesf.path().c_str());
  Heta = dynamic_cast<TH1D*>(fjes_->Get("eta"));
  for(int i = 0 ; i < Heta->GetNbinsX() ; ++i) {
    stringstream hn;
    hn << "down_" << i;
    HptsP.push_back(dynamic_cast<TH1D*>(fjes_->Get(hn.str().c_str())));
  }	

  for(int i = 0 ; i < Heta->GetNbinsX() ; ++i) {
    stringstream hn;
    hn << "up_" << i;
    HptsM.push_back(dynamic_cast<TH1D*>(fjes_->Get(hn.str().c_str())));
  }	

  //Get JER
  fjer_ = new TFile(jerf.path().c_str());
  res_ = dynamic_cast<TF1*>  (fjer_->Get("resolution"));
  sfc_ = dynamic_cast<TH1F*> (fjer_->Get("sf_central"));
  sfu_ = dynamic_cast<TH1F*> (fjer_->Get("sf_up")     );
  sfd_ = dynamic_cast<TH1F*> (fjer_->Get("sf_down")   );

  olddir->cd();    
}
