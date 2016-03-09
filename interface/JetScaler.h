#ifndef HJETSCALER
#define HJETSCALER
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include <TH1D.h>
#include "TF1.h"
#include "TH1F.h"
#include <TDirectory.h>

#include "IDJet.h"
#include "TMath.h"

using namespace std;

class JetScaler {
private:
  TH1D* Heta;
  vector<TH1D*> HptsP;
  vector<TH1D*> HptsM;
  TFile *fjes_, *fjer_;
  
  TF1 *res_;
  TH1F *sfc_;
  TH1F *sfu_;
  TH1F *sfd_;

public:
  static URStreamer* streamer;
  JetScaler(const string filename, const string resfile)
  {
    Init(filename, resfile);
  }

  JetScaler():
    Heta(0),
    HptsP(),
    HptsM(),
    sfc_(),
    sfu_(),
    sfd_() {}

  ~JetScaler() {
    delete fjer_;
    delete fjes_;
  }

  void Init(const string filename, const string resfile) {//like the constructior, but outside it
    TDirectory* olddir = gDirectory;
    fjes_ = new TFile(filename.c_str());
    Heta = dynamic_cast<TH1D*>(fjes_->Get("eta"));
    for(int i = 0 ; i < Heta->GetNbinsX() ; ++i)
    {
      stringstream hn;
      hn << "down_" << i;
      HptsP.push_back(dynamic_cast<TH1D*>(fjes_->Get(hn.str().c_str())));
    }	
    for(int i = 0 ; i < Heta->GetNbinsX() ; ++i)
    {
      stringstream hn;
      hn << "up_" << i;
      HptsM.push_back(dynamic_cast<TH1D*>(fjes_->Get(hn.str().c_str())));
    }	

    fjer_ = new TFile(resfile.c_str());
    res_ = dynamic_cast<TF1*>  (fjer_->Get("resolution"));
    sfc_ = dynamic_cast<TH1F*> (fjer_->Get("sf_central"));
    sfu_ = dynamic_cast<TH1F*> (fjer_->Get("sf_up")     );
    sfd_ = dynamic_cast<TH1F*> (fjer_->Get("sf_down")   );

    olddir->cd();    
  }
  double GetUncP(const IDJet& jet)
  {
    int etabin = Heta->FindFixBin(jet.Eta()) -1;
    int ptbin = HptsP[etabin]->FindFixBin(jet.Pt());
    return(HptsP[etabin]->GetBinContent(ptbin));
  }
  double GetUncM(const IDJet& jet)
  {
    int etabin = Heta->FindFixBin(jet.Eta()) -1;
    int ptbin = HptsM[etabin]->FindFixBin(jet.Pt());
    return(HptsM[etabin]->GetBinContent(ptbin));
  }

  double JERSigmaUp(const IDJet& jet) {
    return jer_sigma(jet, sfu_);
  }

  double JERSigmaDw(const IDJet& jet) {
    return jer_sigma(jet, sfd_);
  }

  double JERSigma(const IDJet& jet) {
    return jer_sigma(jet, sfc_);
  }

private:
  double jer_sigma(const IDJet& jet, TH1F* hsf) {
    double smc = res_->Eval(jet.Pt());
    double sf = hsf->GetBinContent( hsf->GetXaxis()->FindFixBin(jet.Eta()) );
    return TMath::Sqrt(sf*sf-1)*smc;
  }
  // void scale_up(IDJet& jet) {
  //   double sf = 1+GetUncP(jet);
  //   jet.SetPxPyPzE(jet.Px()*sf, jet.Py()*sf, jet.Pz()*sf, jet.E()*sf);
  // }
  // void scale_down(IDJet& jet) {
  //   double sf = 1-GetUncM(jet);
  //   jet.SetPxPyPzE(jet.Px()*sf, jet.Py()*sf, jet.Pz()*sf, jet.E()*sf);
  // }
};

#endif
