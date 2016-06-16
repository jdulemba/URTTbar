#ifndef HJETSCALER
#define HJETSCALER
#include <vector>

#include <TH1D.h>
#include "TF1.h"
#include "TH1F.h"
#include <TDirectory.h>

#include "Analyses/URTTbar/interface/IDJet.h"
#include "TMath.h"
#include <memory>

class JetScaler {
private:
  std::shared_ptr<TH1D> Heta;
  std::vector<std::shared_ptr<TH1D> > HptsP;
  std::vector<std::shared_ptr<TH1D> > HptsM;
  //TFile *fjes_, *fjer_;

  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
  std::shared_ptr<TF1 > res_;
  std::shared_ptr<TH1F> sfc_;
  std::shared_ptr<TH1F> sfu_;
  std::shared_ptr<TH1F> sfd_;

  template <class T>
  std::shared_ptr<T> get_from(TFile &file, std::string path, std::string newname) {
    T* original = (T*) file.Get( path.c_str() );
    if(!original) {
      Logger::log().fatal() << "Could not get " << path << " from the file " << file.GetName() << "!" << std::endl;
      throw 42;
    }
    std::shared_ptr<T> ptr((T*) original->Clone(newname.c_str()));
    return ptr;
  }

public:
  static JetScaler& instance() {
    static JetScaler val;
    return val; 
  }

  ~JetScaler() {
    // delete fjer_;
    // delete fjes_;
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

  inline double jer_sf(const IDJet& jet, std::shared_ptr<TH1F> hsf) {
    //allows easier refactoring when JER changes
    return hsf->GetBinContent( hsf->GetXaxis()->FindFixBin(jet.Eta()) );
  }

  double jer_sigma(const IDJet& jet, std::shared_ptr<TH1F> hsf) {
    double smc = mc_resolution(jet);
    double sf  = jer_sf(jet, hsf);
    return TMath::Sqrt(sf*sf-1)*smc;
  }
};

#endif
