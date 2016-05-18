#ifndef HJETSCALER
#define HJETSCALER
#include <vector>

#include <TH1D.h>
#include "TF1.h"
#include "TH1F.h"
#include <TDirectory.h>

#include "Analyses/URTTbar/interface/IDJet.h"
#include "TMath.h"

class JetScaler {
private:
  TH1D* Heta;
  std::vector<TH1D*> HptsP;
  std::vector<TH1D*> HptsM;
  TFile *fjes_, *fjer_;

  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
  TF1 *res_;
  TH1F *sfc_;
  TH1F *sfu_;
  TH1F *sfd_;

public:
  static JetScaler& instance() {
    static JetScaler val;
    return val; 
  }

  ~JetScaler() {
    delete fjer_;
    delete fjes_;
  }

  double GetUncP(const IDJet& jet) {
    int etabin = Heta->FindFixBin(jet.Eta()) -1;
    int ptbin = HptsP[etabin]->FindFixBin(jet.Pt());
    return(HptsP[etabin]->GetBinContent(ptbin));
  }
  double GetUncM(const IDJet& jet) {
    int etabin = Heta->FindFixBin(jet.Eta()) -1;
    int ptbin = HptsM[etabin]->FindFixBin(jet.Pt());
    return(HptsM[etabin]->GetBinContent(ptbin));
  }

  //returns absolute resolution
  double resolution(const IDJet& jet) {return jet.E()*mc_resolution(jet)*jer_sf(jet, sfc_);}

  double JERSigmaUp(const IDJet& jet) {return jer_sigma(jet, sfu_);}
  double JERSigmaDw(const IDJet& jet) {return jer_sigma(jet, sfd_);}
  double JERSigma(const IDJet& jet) {return jer_sigma(jet, sfc_);}

private:
  JetScaler();  

  //No need to implement those
  JetScaler(JetScaler const&);
  void operator=(JetScaler const&);

  inline double mc_resolution(const IDJet& jet) {
    //allows easier refactoring when JER changes
    return res_->Eval(jet.Pt());
  }

  inline double jer_sf(const IDJet& jet, TH1F* hsf) {
    //allows easier refactoring when JER changes
    return hsf->GetBinContent( hsf->GetXaxis()->FindFixBin(jet.Eta()) );
  }

  double jer_sigma(const IDJet& jet, TH1F* hsf) {
    double smc = mc_resolution(jet);
    double sf  = jer_sf(jet, hsf);
    return TMath::Sqrt(sf*sf-1)*smc;
  }
};

#endif
