#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/Wards.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include "TH1I.h"

using namespace TMath;
LeptonSF::LeptonSF(std::string parname, bool ptx):
  pt_as_x_(ptx)
{
  URParser &parser = URParser::instance();
  parser.addCfgParameter<std::string>("general", parname, "lepton scale factor file");
  
  parser.parseArguments();

  DataFile sf_file = parser.getCfgPar<std::string>("general", parname);
  TFile file(sf_file.path().c_str());
	HistoOwnershipWard hward;
	DirectoryWard dward;
  id_   = get_from<TH2D>(file, "id" , "id_clone");
  iso_  = get_from<TH2D>(file, "iso", "iso_clone");
  trig_ = get_from<TH2D>(file, "trg", "trg_clone");
	trk_  = get_from<TH2D>(file, "trk", "trk_clone");
  
  TH1I *h = (TH1I*) file.Get("info");
  if(h) {
    pt_as_x_ = (h->GetBinContent(0) > 0.5);
    for(int i=0; i<4; i++)
      abs_etas_[i] = (h->GetBinContent(i+1) > 0.5);
  }
}

double LeptonSF::get_2d_weight(std::shared_ptr<TH2> h, double pt, double eta) const {
  if(!h.get()) return 1.;
  int pt_bin = h->GetXaxis()->FindFixBin((pt_as_x_) ? pt : eta);
  pt_bin = TMath::Max(pt_bin, 1);
  pt_bin = TMath::Min(pt_bin, h->GetNbinsX());

  int eta_bin = h->GetYaxis()->FindFixBin((pt_as_x_) ? eta : pt);
  eta_bin = TMath::Max(eta_bin, 1);
  eta_bin = TMath::Min(eta_bin, h->GetNbinsY());

  return h->GetBinContent(pt_bin, eta_bin);
}

double LeptonSF::get_1d_weight(std::shared_ptr<TH1> h, double pt, double eta) const {
  if(!h.get()) return 1.;
  int pt_bin = h->GetXaxis()->FindFixBin((pt_as_x_) ? pt : eta);
  pt_bin = TMath::Max(pt_bin, 1);
  pt_bin = TMath::Min(pt_bin, h->GetNbinsX());

  return h->GetBinContent(pt_bin);
}

double LeptonSF::get_sf(double pt, double eta) const {
  return  get_2d_weight(trig_, pt, abs_etas_[0] ? Abs(eta) : eta) *
    get_2d_weight(id_ , pt, abs_etas_[1] ? Abs(eta) : eta) *
    get_2d_weight(iso_, pt, abs_etas_[2] ? Abs(eta) : eta) *
		get_1d_weight(trk_, pt, abs_etas_[3] ? Abs(eta) : eta);
}
