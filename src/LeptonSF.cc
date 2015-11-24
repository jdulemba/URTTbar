#include "LeptonSF.h"
#include "DataFile.h"
#include "URParser.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>

LeptonSF::LeptonSF(std::string parname, bool ptx):
  pt_as_x_(ptx)
{
  URParser &parser = URParser::instance();
  parser.addCfgParameter<std::string>("general", parname, "lepton scale factor file");
  
  parser.parseArguments();

  DataFile sf_file = parser.getCfgPar<std::string>("general", parname);
  TFile file(sf_file.path().c_str());
  TH1::AddDirectory(false); //TODO: add ward and raise warning if histo is missing
  TH2F *ptr = (TH2F*) file.Get("id");
  id_   = (ptr) ? (TH2F*) ptr->Clone("id_clone") : NULL;
  ptr = (TH2F*) file.Get("iso");
  iso_  = (ptr) ? (TH2F*) ptr->Clone("iso_clone") : NULL;
  ptr = (TH2F*) file.Get("trg");
  trig_ = (ptr) ? (TH2F*) ptr->Clone("trg_clone") : NULL;
  TH1::AddDirectory(true);
}

double LeptonSF::get_weight(TH2 *h, double pt, double eta) const {
  if(!h) return 1.;
  int pt_bin = h->GetXaxis()->FindFixBin((pt_as_x_) ? pt : eta);
  pt_bin = TMath::Max(pt_bin, 1);
  pt_bin = TMath::Min(pt_bin, h->GetNbinsX());

  int eta_bin = h->GetYaxis()->FindFixBin((pt_as_x_) ? eta : pt);
  eta_bin = TMath::Max(eta_bin, 1);
  eta_bin = TMath::Min(eta_bin, h->GetNbinsY());

  return h->GetBinContent(pt_bin, eta_bin);
}

double LeptonSF::get_sf(double pt, double eta) const {
  return get_weight(id_, pt, eta) *
    get_weight(iso_, pt, eta) *
    get_weight(trig_, pt, eta);
}
