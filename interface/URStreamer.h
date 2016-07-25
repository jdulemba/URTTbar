#ifndef URStreamer_h
#define URStreamer_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include <vector>
using namespace std;

class Npnlolhe{
friend class URStreamer;
public:
//  Npnlolhe(const Int_t &i_npnlo_):
//    
//  {}
  Npnlolhe():
    npnlo_(0)
  {}
  Int_t npnlo() const {return npnlo_;}
private:
  Int_t npnlo_;
  void setnpnlo(const Int_t value) {npnlo_ = value;}
};

class Met{
friend class URStreamer;
public:
//  Met(const float &i_px_,const float &i_py_,const float &i_pxsmear_,const float &i_pysmear_,const float &i_pxunc_,const float &i_pyunc_,const float &i_pxuncJES_,const float &i_pyuncJES_,const float &i_pxuncJER_,const float &i_pyuncJER_):
//    
//  {}
  Met():
    px_(0),
    py_(0),
    pxsmear_(0),
    pysmear_(0),
    pxunc_(0),
    pyunc_(0),
    pxuncJES_(0),
    pyuncJES_(0),
    pxuncJER_(0),
    pyuncJER_(0)
  {}
  float px() const {return px_;}
  float py() const {return py_;}
  float pxsmear() const {return pxsmear_;}
  float pysmear() const {return pysmear_;}
  float pxunc() const {return pxunc_;}
  float pyunc() const {return pyunc_;}
  float pxuncJES() const {return pxuncJES_;}
  float pyuncJES() const {return pyuncJES_;}
  float pxuncJER() const {return pxuncJER_;}
  float pyuncJER() const {return pyuncJER_;}
private:
  float px_;
  float py_;
  float pxsmear_;
  float pysmear_;
  float pxunc_;
  float pyunc_;
  float pxuncJES_;
  float pyuncJES_;
  float pxuncJER_;
  float pyuncJER_;
  void setpx(const float value) {px_ = value;}
  void setpy(const float value) {py_ = value;}
  void setpxsmear(const float value) {pxsmear_ = value;}
  void setpysmear(const float value) {pysmear_ = value;}
  void setpxunc(const float value) {pxunc_ = value;}
  void setpyunc(const float value) {pyunc_ = value;}
  void setpxuncJES(const float value) {pxuncJES_ = value;}
  void setpyuncJES(const float value) {pyuncJES_ = value;}
  void setpxuncJER(const float value) {pxuncJER_ = value;}
  void setpyuncJER(const float value) {pyuncJER_ = value;}
};

class Mcweight{
friend class URStreamer;
public:
//  Mcweight(const float &i_weights_):
//    
//  {}
  Mcweight():
    weights_(0)
  {}
  float weights() const {return weights_;}
private:
  float weights_;
  void setweights(const float value) {weights_ = value;}
};

class Geninfo{
friend class URStreamer;
public:
//  Geninfo(const Float_t &i_weight_,const Float_t &i_pdfid1_,const Float_t &i_pdfid2_,const Float_t &i_x1_,const Float_t &i_x2_,const Float_t &i_renScale_):
//    
//  {}
  Geninfo():
    weight_(0),
    pdfid1_(0),
    pdfid2_(0),
    x1_(0),
    x2_(0),
    renScale_(0)
  {}
  Float_t weight() const {return weight_;}
  Float_t pdfid1() const {return pdfid1_;}
  Float_t pdfid2() const {return pdfid2_;}
  Float_t x1() const {return x1_;}
  Float_t x2() const {return x2_;}
  Float_t renScale() const {return renScale_;}
private:
  Float_t weight_;
  Float_t pdfid1_;
  Float_t pdfid2_;
  Float_t x1_;
  Float_t x2_;
  Float_t renScale_;
  void setweight(const Float_t value) {weight_ = value;}
  void setpdfid1(const Float_t value) {pdfid1_ = value;}
  void setpdfid2(const Float_t value) {pdfid2_ = value;}
  void setx1(const Float_t value) {x1_ = value;}
  void setx2(const Float_t value) {x2_ = value;}
  void setrenScale(const Float_t value) {renScale_ = value;}
};

class Photon: public TLorentzVector{
friend class URStreamer;
public:
//  Photon(const int &i_charge_,const float &i_x_,const float &i_y_,const float &i_z_,const float &i_energy_,const float &i_rawEnergy_,const float &i_phiWidth_,const float &i_etaWidth_,const float &i_e3x3_,const float &i_maxCrystalEnergy_,const bool &i_isEB_,const bool &i_isEE_,const bool &i_isPFlowPhoton_,const bool &i_hasConversionTracks_,const bool &i_hasPixelSeed_):
//    
//  {}
  Photon():
    TLorentzVector(),
    charge_(0),
    x_(0),
    y_(0),
    z_(0),
    energy_(0),
    rawEnergy_(0),
    phiWidth_(0),
    etaWidth_(0),
    e3x3_(0),
    maxCrystalEnergy_(0),
    isEB_(0),
    isEE_(0),
    isPFlowPhoton_(0),
    hasConversionTracks_(0),
    hasPixelSeed_(0)
  {}
  int charge() const {return charge_;}
  float x() const {return x_;}
  float y() const {return y_;}
  float z() const {return z_;}
  float energy() const {return energy_;}
  float rawEnergy() const {return rawEnergy_;}
  float phiWidth() const {return phiWidth_;}
  float etaWidth() const {return etaWidth_;}
  float e3x3() const {return e3x3_;}
  float maxCrystalEnergy() const {return maxCrystalEnergy_;}
  bool isEB() const {return isEB_;}
  bool isEE() const {return isEE_;}
  bool isPFlowPhoton() const {return isPFlowPhoton_;}
  bool hasConversionTracks() const {return hasConversionTracks_;}
  bool hasPixelSeed() const {return hasPixelSeed_;}
private:
  int charge_;
  float x_;
  float y_;
  float z_;
  float energy_;
  float rawEnergy_;
  float phiWidth_;
  float etaWidth_;
  float e3x3_;
  float maxCrystalEnergy_;
  bool isEB_;
  bool isEE_;
  bool isPFlowPhoton_;
  bool hasConversionTracks_;
  bool hasPixelSeed_;
  void setcharge(const int value) {charge_ = value;}
  void setx(const float value) {x_ = value;}
  void sety(const float value) {y_ = value;}
  void setz(const float value) {z_ = value;}
  void setenergy(const float value) {energy_ = value;}
  void setrawEnergy(const float value) {rawEnergy_ = value;}
  void setphiWidth(const float value) {phiWidth_ = value;}
  void setetaWidth(const float value) {etaWidth_ = value;}
  void sete3x3(const float value) {e3x3_ = value;}
  void setmaxCrystalEnergy(const float value) {maxCrystalEnergy_ = value;}
  void setisEB(const bool value) {isEB_ = value;}
  void setisEE(const bool value) {isEE_ = value;}
  void setisPFlowPhoton(const bool value) {isPFlowPhoton_ = value;}
  void sethasConversionTracks(const bool value) {hasConversionTracks_ = value;}
  void sethasPixelSeed(const bool value) {hasPixelSeed_ = value;}
  void setLotentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Lhepaticle{
friend class URStreamer;
public:
//  Lhepaticle(const float &i_px_,const float &i_py_,const float &i_pz_,const float &i_e_,const int &i_pdgid_,const int &i_status_,const int &i_fmother_,const int &i_lmother_):
//    
//  {}
  Lhepaticle():
    px_(0),
    py_(0),
    pz_(0),
    e_(0),
    pdgid_(0),
    status_(0),
    fmother_(0),
    lmother_(0)
  {}
  float px() const {return px_;}
  float py() const {return py_;}
  float pz() const {return pz_;}
  float e() const {return e_;}
  int pdgid() const {return pdgid_;}
  int status() const {return status_;}
  int fmother() const {return fmother_;}
  int lmother() const {return lmother_;}
private:
  float px_;
  float py_;
  float pz_;
  float e_;
  int pdgid_;
  int status_;
  int fmother_;
  int lmother_;
  void setpx(const float value) {px_ = value;}
  void setpy(const float value) {py_ = value;}
  void setpz(const float value) {pz_ = value;}
  void sete(const float value) {e_ = value;}
  void setpdgid(const int value) {pdgid_ = value;}
  void setstatus(const int value) {status_ = value;}
  void setfmother(const int value) {fmother_ = value;}
  void setlmother(const int value) {lmother_ = value;}
};

class Filter{
friend class URStreamer;
public:
//  Filter(const Int_t &i_Flag_goodVertices_,const Int_t &i_Flag_CSCTightHaloFilter_,const Int_t &i_Flag_trkPOGFilters_,const Int_t &i_Flag_trkPOG_logErrorTooManyClusters_,const Int_t &i_Flag_EcalDeadCellTriggerPrimitiveFilter_,const Int_t &i_Flag_ecalLaserCorrFilter_,const Int_t &i_Flag_trkPOG_manystripclus53X_,const Int_t &i_Flag_eeBadScFilter_,const Int_t &i_Flag_METFilters_,const Int_t &i_Flag_HBHENoiseFilter_,const Int_t &i_Flag_trkPOG_toomanystripclus53X_,const Int_t &i_Flag_hcalLaserEventFilter_,const Int_t &i_Flag_HBHENoiseIsoFilter_,const Int_t &i_Flag_CSCTightHalo2015Filter_):
//    
//  {}
  Filter():
    Flag_goodVertices_(0),
    Flag_CSCTightHaloFilter_(0),
    Flag_trkPOGFilters_(0),
    Flag_trkPOG_logErrorTooManyClusters_(0),
    Flag_EcalDeadCellTriggerPrimitiveFilter_(0),
    Flag_ecalLaserCorrFilter_(0),
    Flag_trkPOG_manystripclus53X_(0),
    Flag_eeBadScFilter_(0),
    Flag_METFilters_(0),
    Flag_HBHENoiseFilter_(0),
    Flag_trkPOG_toomanystripclus53X_(0),
    Flag_hcalLaserEventFilter_(0),
    Flag_HBHENoiseIsoFilter_(0),
    Flag_CSCTightHalo2015Filter_(0)
  {}
  Int_t Flag_goodVertices() const {return Flag_goodVertices_;}
  Int_t Flag_CSCTightHaloFilter() const {return Flag_CSCTightHaloFilter_;}
  Int_t Flag_trkPOGFilters() const {return Flag_trkPOGFilters_;}
  Int_t Flag_trkPOG_logErrorTooManyClusters() const {return Flag_trkPOG_logErrorTooManyClusters_;}
  Int_t Flag_EcalDeadCellTriggerPrimitiveFilter() const {return Flag_EcalDeadCellTriggerPrimitiveFilter_;}
  Int_t Flag_ecalLaserCorrFilter() const {return Flag_ecalLaserCorrFilter_;}
  Int_t Flag_trkPOG_manystripclus53X() const {return Flag_trkPOG_manystripclus53X_;}
  Int_t Flag_eeBadScFilter() const {return Flag_eeBadScFilter_;}
  Int_t Flag_METFilters() const {return Flag_METFilters_;}
  Int_t Flag_HBHENoiseFilter() const {return Flag_HBHENoiseFilter_;}
  Int_t Flag_trkPOG_toomanystripclus53X() const {return Flag_trkPOG_toomanystripclus53X_;}
  Int_t Flag_hcalLaserEventFilter() const {return Flag_hcalLaserEventFilter_;}
  Int_t Flag_HBHENoiseIsoFilter() const {return Flag_HBHENoiseIsoFilter_;}
  Int_t Flag_CSCTightHalo2015Filter() const {return Flag_CSCTightHalo2015Filter_;}
private:
  Int_t Flag_goodVertices_;
  Int_t Flag_CSCTightHaloFilter_;
  Int_t Flag_trkPOGFilters_;
  Int_t Flag_trkPOG_logErrorTooManyClusters_;
  Int_t Flag_EcalDeadCellTriggerPrimitiveFilter_;
  Int_t Flag_ecalLaserCorrFilter_;
  Int_t Flag_trkPOG_manystripclus53X_;
  Int_t Flag_eeBadScFilter_;
  Int_t Flag_METFilters_;
  Int_t Flag_HBHENoiseFilter_;
  Int_t Flag_trkPOG_toomanystripclus53X_;
  Int_t Flag_hcalLaserEventFilter_;
  Int_t Flag_HBHENoiseIsoFilter_;
  Int_t Flag_CSCTightHalo2015Filter_;
  void setFlag_goodVertices(const Int_t value) {Flag_goodVertices_ = value;}
  void setFlag_CSCTightHaloFilter(const Int_t value) {Flag_CSCTightHaloFilter_ = value;}
  void setFlag_trkPOGFilters(const Int_t value) {Flag_trkPOGFilters_ = value;}
  void setFlag_trkPOG_logErrorTooManyClusters(const Int_t value) {Flag_trkPOG_logErrorTooManyClusters_ = value;}
  void setFlag_EcalDeadCellTriggerPrimitiveFilter(const Int_t value) {Flag_EcalDeadCellTriggerPrimitiveFilter_ = value;}
  void setFlag_ecalLaserCorrFilter(const Int_t value) {Flag_ecalLaserCorrFilter_ = value;}
  void setFlag_trkPOG_manystripclus53X(const Int_t value) {Flag_trkPOG_manystripclus53X_ = value;}
  void setFlag_eeBadScFilter(const Int_t value) {Flag_eeBadScFilter_ = value;}
  void setFlag_METFilters(const Int_t value) {Flag_METFilters_ = value;}
  void setFlag_HBHENoiseFilter(const Int_t value) {Flag_HBHENoiseFilter_ = value;}
  void setFlag_trkPOG_toomanystripclus53X(const Int_t value) {Flag_trkPOG_toomanystripclus53X_ = value;}
  void setFlag_hcalLaserEventFilter(const Int_t value) {Flag_hcalLaserEventFilter_ = value;}
  void setFlag_HBHENoiseIsoFilter(const Int_t value) {Flag_HBHENoiseIsoFilter_ = value;}
  void setFlag_CSCTightHalo2015Filter(const Int_t value) {Flag_CSCTightHalo2015Filter_ = value;}
};

class Trigger{
friend class URStreamer;
public:
//  Trigger(const Int_t &i_HLT_Ele27_eta2p1_WPLoose_Gsf_,const Int_t &i_HLT_DoubleIsoMu17_eta2p1_,const Int_t &i_HLT_notexists_,const Int_t &i_HLT_IsoMu24_eta2p1_,const Int_t &i_HLT_IsoMu20_eta2p1_,const Int_t &i_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_,const Int_t &i_HLT_IsoMu18_,const Int_t &i_HLT_IsoTkMu20_,const Int_t &i_HLT_IsoMu20_,const Int_t &i_HLT_IsoTkMu27_,const Int_t &i_HLT_IsoMu27_,const Int_t &i_HLT_Mu50_,const Int_t &i_HLT_Mu45_eta2p1_,const Int_t &i_HLT_Ele22_eta2p1_WPLoose_Gsf_,const Int_t &i_HLT_Ele23_WPLoose_Gsf_,const Int_t &i_HLT_Ele27_WPLoose_Gsf_,const Int_t &i_HLT_IsoMu22_,const Int_t &i_HLT_IsoTkMu22_):
//    
//  {}
  Trigger():
    HLT_Ele27_eta2p1_WPLoose_Gsf_(0),
    HLT_DoubleIsoMu17_eta2p1_(0),
    HLT_notexists_(0),
    HLT_IsoMu24_eta2p1_(0),
    HLT_IsoMu20_eta2p1_(0),
    HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_(0),
    HLT_IsoMu18_(0),
    HLT_IsoTkMu20_(0),
    HLT_IsoMu20_(0),
    HLT_IsoTkMu27_(0),
    HLT_IsoMu27_(0),
    HLT_Mu50_(0),
    HLT_Mu45_eta2p1_(0),
    HLT_Ele22_eta2p1_WPLoose_Gsf_(0),
    HLT_Ele23_WPLoose_Gsf_(0),
    HLT_Ele27_WPLoose_Gsf_(0),
    HLT_IsoMu22_(0),
    HLT_IsoTkMu22_(0)
  {}
  Int_t HLT_Ele27_eta2p1_WPLoose_Gsf() const {return HLT_Ele27_eta2p1_WPLoose_Gsf_;}
  Int_t HLT_DoubleIsoMu17_eta2p1() const {return HLT_DoubleIsoMu17_eta2p1_;}
  Int_t HLT_notexists() const {return HLT_notexists_;}
  Int_t HLT_IsoMu24_eta2p1() const {return HLT_IsoMu24_eta2p1_;}
  Int_t HLT_IsoMu20_eta2p1() const {return HLT_IsoMu20_eta2p1_;}
  Int_t HLT_DoubleEle33_CaloIdL_GsfTrkIdVL() const {return HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_;}
  Int_t HLT_IsoMu18() const {return HLT_IsoMu18_;}
  Int_t HLT_IsoTkMu20() const {return HLT_IsoTkMu20_;}
  Int_t HLT_IsoMu20() const {return HLT_IsoMu20_;}
  Int_t HLT_IsoTkMu27() const {return HLT_IsoTkMu27_;}
  Int_t HLT_IsoMu27() const {return HLT_IsoMu27_;}
  Int_t HLT_Mu50() const {return HLT_Mu50_;}
  Int_t HLT_Mu45_eta2p1() const {return HLT_Mu45_eta2p1_;}
  Int_t HLT_Ele22_eta2p1_WPLoose_Gsf() const {return HLT_Ele22_eta2p1_WPLoose_Gsf_;}
  Int_t HLT_Ele23_WPLoose_Gsf() const {return HLT_Ele23_WPLoose_Gsf_;}
  Int_t HLT_Ele27_WPLoose_Gsf() const {return HLT_Ele27_WPLoose_Gsf_;}
  Int_t HLT_IsoMu22() const {return HLT_IsoMu22_;}
  Int_t HLT_IsoTkMu22() const {return HLT_IsoTkMu22_;}
private:
  Int_t HLT_Ele27_eta2p1_WPLoose_Gsf_;
  Int_t HLT_DoubleIsoMu17_eta2p1_;
  Int_t HLT_notexists_;
  Int_t HLT_IsoMu24_eta2p1_;
  Int_t HLT_IsoMu20_eta2p1_;
  Int_t HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_;
  Int_t HLT_IsoMu18_;
  Int_t HLT_IsoTkMu20_;
  Int_t HLT_IsoMu20_;
  Int_t HLT_IsoTkMu27_;
  Int_t HLT_IsoMu27_;
  Int_t HLT_Mu50_;
  Int_t HLT_Mu45_eta2p1_;
  Int_t HLT_Ele22_eta2p1_WPLoose_Gsf_;
  Int_t HLT_Ele23_WPLoose_Gsf_;
  Int_t HLT_Ele27_WPLoose_Gsf_;
  Int_t HLT_IsoMu22_;
  Int_t HLT_IsoTkMu22_;
  void setHLT_Ele27_eta2p1_WPLoose_Gsf(const Int_t value) {HLT_Ele27_eta2p1_WPLoose_Gsf_ = value;}
  void setHLT_DoubleIsoMu17_eta2p1(const Int_t value) {HLT_DoubleIsoMu17_eta2p1_ = value;}
  void setHLT_notexists(const Int_t value) {HLT_notexists_ = value;}
  void setHLT_IsoMu24_eta2p1(const Int_t value) {HLT_IsoMu24_eta2p1_ = value;}
  void setHLT_IsoMu20_eta2p1(const Int_t value) {HLT_IsoMu20_eta2p1_ = value;}
  void setHLT_DoubleEle33_CaloIdL_GsfTrkIdVL(const Int_t value) {HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_ = value;}
  void setHLT_IsoMu18(const Int_t value) {HLT_IsoMu18_ = value;}
  void setHLT_IsoTkMu20(const Int_t value) {HLT_IsoTkMu20_ = value;}
  void setHLT_IsoMu20(const Int_t value) {HLT_IsoMu20_ = value;}
  void setHLT_IsoTkMu27(const Int_t value) {HLT_IsoTkMu27_ = value;}
  void setHLT_IsoMu27(const Int_t value) {HLT_IsoMu27_ = value;}
  void setHLT_Mu50(const Int_t value) {HLT_Mu50_ = value;}
  void setHLT_Mu45_eta2p1(const Int_t value) {HLT_Mu45_eta2p1_ = value;}
  void setHLT_Ele22_eta2p1_WPLoose_Gsf(const Int_t value) {HLT_Ele22_eta2p1_WPLoose_Gsf_ = value;}
  void setHLT_Ele23_WPLoose_Gsf(const Int_t value) {HLT_Ele23_WPLoose_Gsf_ = value;}
  void setHLT_Ele27_WPLoose_Gsf(const Int_t value) {HLT_Ele27_WPLoose_Gsf_ = value;}
  void setHLT_IsoMu22(const Int_t value) {HLT_IsoMu22_ = value;}
  void setHLT_IsoTkMu22(const Int_t value) {HLT_IsoTkMu22_ = value;}
};

class Electron: public TLorentzVector{
friend class URStreamer;
public:
//  Electron(const int &i_charge_,const float &i_chargedIso_,const float &i_neutralIso_,const float &i_photonIso_,const float &i_puIso_,const float &i_dB_,const float &i_ipDXY_,const float &i_dz_,const float &i_nMissingInnerHits_,const float &i_r9_,const float &i_ESCOverETrack_,const float &i_DEtaSCTrk_,const float &i_DPhiSCTrk_,const float &i_ecalEnergy_,const bool &i_passConversionVeto_,const bool &i_isEB_,const bool &i_isEE_,const bool &i_isEBGap_,const bool &i_isEBEtaGap_,const bool &i_isEBPhiGap_,const bool &i_isEEGap_,const bool &i_isEERingGap_,const bool &i_isEEDeeGap_,const bool &i_isEBEEGap_,const bool &i_isElectron_,const bool &i_ecalSeed_,const bool &i_trackSeed_,const float &i_eidCutLoose_,const float &i_eidCutMedium_,const float &i_eidCutTight_,const float &i_eidCutVeto_,const float &i_eidMVAWP80_,const float &i_eidMVAWP90_,const float &i_eidTrgMVAWP80_,const float &i_eidTrgMVAWP90_,const float &i_pfHadronIso_,const float &i_pfNeutralIso_,const float &i_pfPhotonIso_,const bool &i_HLT_Ele27_eta2p1_WPLoose_Gsf_,const bool &i_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_,const bool &i_HLT_Ele22_eta2p1_WPLoose_Gsf_,const bool &i_HLT_Ele23_WPLoose_Gsf_,const bool &i_HLT_Ele27_WPLoose_Gsf_,const float &i_e1x5_,const float &i_e5x5_,const float &i_sigmaIEtaIEta_,const float &i_full5x5_sigmaIEtaIEta_,const float &i_sigmaIPhiIPhi_,const float &i_hadronicOverEM_,const float &i_x_,const float &i_y_,const float &i_z_,const float &i_energy_,const float &i_rawEnergy_,const float &i_phiWidth_,const float &i_etaWidth_):
//    
//  {}
  Electron():
    TLorentzVector(),
    charge_(0),
    chargedIso_(0),
    neutralIso_(0),
    photonIso_(0),
    puIso_(0),
    dB_(0),
    ipDXY_(0),
    dz_(0),
    nMissingInnerHits_(0),
    r9_(0),
    ESCOverETrack_(0),
    DEtaSCTrk_(0),
    DPhiSCTrk_(0),
    ecalEnergy_(0),
    passConversionVeto_(0),
    isEB_(0),
    isEE_(0),
    isEBGap_(0),
    isEBEtaGap_(0),
    isEBPhiGap_(0),
    isEEGap_(0),
    isEERingGap_(0),
    isEEDeeGap_(0),
    isEBEEGap_(0),
    isElectron_(0),
    ecalSeed_(0),
    trackSeed_(0),
    eidCutLoose_(0),
    eidCutMedium_(0),
    eidCutTight_(0),
    eidCutVeto_(0),
    eidMVAWP80_(0),
    eidMVAWP90_(0),
    eidTrgMVAWP80_(0),
    eidTrgMVAWP90_(0),
    pfHadronIso_(0),
    pfNeutralIso_(0),
    pfPhotonIso_(0),
    HLT_Ele27_eta2p1_WPLoose_Gsf_(0),
    HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_(0),
    HLT_Ele22_eta2p1_WPLoose_Gsf_(0),
    HLT_Ele23_WPLoose_Gsf_(0),
    HLT_Ele27_WPLoose_Gsf_(0),
    e1x5_(0),
    e5x5_(0),
    sigmaIEtaIEta_(0),
    full5x5_sigmaIEtaIEta_(0),
    sigmaIPhiIPhi_(0),
    hadronicOverEM_(0),
    x_(0),
    y_(0),
    z_(0),
    energy_(0),
    rawEnergy_(0),
    phiWidth_(0),
    etaWidth_(0)
  {}
  int charge() const {return charge_;}
  float chargedIso() const {return chargedIso_;}
  float neutralIso() const {return neutralIso_;}
  float photonIso() const {return photonIso_;}
  float puIso() const {return puIso_;}
  float dB() const {return dB_;}
  float ipDXY() const {return ipDXY_;}
  float dz() const {return dz_;}
  float nMissingInnerHits() const {return nMissingInnerHits_;}
  float r9() const {return r9_;}
  float ESCOverETrack() const {return ESCOverETrack_;}
  float DEtaSCTrk() const {return DEtaSCTrk_;}
  float DPhiSCTrk() const {return DPhiSCTrk_;}
  float ecalEnergy() const {return ecalEnergy_;}
  bool passConversionVeto() const {return passConversionVeto_;}
  bool isEB() const {return isEB_;}
  bool isEE() const {return isEE_;}
  bool isEBGap() const {return isEBGap_;}
  bool isEBEtaGap() const {return isEBEtaGap_;}
  bool isEBPhiGap() const {return isEBPhiGap_;}
  bool isEEGap() const {return isEEGap_;}
  bool isEERingGap() const {return isEERingGap_;}
  bool isEEDeeGap() const {return isEEDeeGap_;}
  bool isEBEEGap() const {return isEBEEGap_;}
  bool isElectron() const {return isElectron_;}
  bool ecalSeed() const {return ecalSeed_;}
  bool trackSeed() const {return trackSeed_;}
  float eidCutLoose() const {return eidCutLoose_;}
  float eidCutMedium() const {return eidCutMedium_;}
  float eidCutTight() const {return eidCutTight_;}
  float eidCutVeto() const {return eidCutVeto_;}
  float eidMVAWP80() const {return eidMVAWP80_;}
  float eidMVAWP90() const {return eidMVAWP90_;}
  float eidTrgMVAWP80() const {return eidTrgMVAWP80_;}
  float eidTrgMVAWP90() const {return eidTrgMVAWP90_;}
  float pfHadronIso() const {return pfHadronIso_;}
  float pfNeutralIso() const {return pfNeutralIso_;}
  float pfPhotonIso() const {return pfPhotonIso_;}
  bool HLT_Ele27_eta2p1_WPLoose_Gsf() const {return HLT_Ele27_eta2p1_WPLoose_Gsf_;}
  bool HLT_DoubleEle33_CaloIdL_GsfTrkIdVL() const {return HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_;}
  bool HLT_Ele22_eta2p1_WPLoose_Gsf() const {return HLT_Ele22_eta2p1_WPLoose_Gsf_;}
  bool HLT_Ele23_WPLoose_Gsf() const {return HLT_Ele23_WPLoose_Gsf_;}
  bool HLT_Ele27_WPLoose_Gsf() const {return HLT_Ele27_WPLoose_Gsf_;}
  float e1x5() const {return e1x5_;}
  float e5x5() const {return e5x5_;}
  float sigmaIEtaIEta() const {return sigmaIEtaIEta_;}
  float full5x5_sigmaIEtaIEta() const {return full5x5_sigmaIEtaIEta_;}
  float sigmaIPhiIPhi() const {return sigmaIPhiIPhi_;}
  float hadronicOverEM() const {return hadronicOverEM_;}
  float x() const {return x_;}
  float y() const {return y_;}
  float z() const {return z_;}
  float energy() const {return energy_;}
  float rawEnergy() const {return rawEnergy_;}
  float phiWidth() const {return phiWidth_;}
  float etaWidth() const {return etaWidth_;}
private:
  int charge_;
  float chargedIso_;
  float neutralIso_;
  float photonIso_;
  float puIso_;
  float dB_;
  float ipDXY_;
  float dz_;
  float nMissingInnerHits_;
  float r9_;
  float ESCOverETrack_;
  float DEtaSCTrk_;
  float DPhiSCTrk_;
  float ecalEnergy_;
  bool passConversionVeto_;
  bool isEB_;
  bool isEE_;
  bool isEBGap_;
  bool isEBEtaGap_;
  bool isEBPhiGap_;
  bool isEEGap_;
  bool isEERingGap_;
  bool isEEDeeGap_;
  bool isEBEEGap_;
  bool isElectron_;
  bool ecalSeed_;
  bool trackSeed_;
  float eidCutLoose_;
  float eidCutMedium_;
  float eidCutTight_;
  float eidCutVeto_;
  float eidMVAWP80_;
  float eidMVAWP90_;
  float eidTrgMVAWP80_;
  float eidTrgMVAWP90_;
  float pfHadronIso_;
  float pfNeutralIso_;
  float pfPhotonIso_;
  bool HLT_Ele27_eta2p1_WPLoose_Gsf_;
  bool HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_;
  bool HLT_Ele22_eta2p1_WPLoose_Gsf_;
  bool HLT_Ele23_WPLoose_Gsf_;
  bool HLT_Ele27_WPLoose_Gsf_;
  float e1x5_;
  float e5x5_;
  float sigmaIEtaIEta_;
  float full5x5_sigmaIEtaIEta_;
  float sigmaIPhiIPhi_;
  float hadronicOverEM_;
  float x_;
  float y_;
  float z_;
  float energy_;
  float rawEnergy_;
  float phiWidth_;
  float etaWidth_;
  void setcharge(const int value) {charge_ = value;}
  void setchargedIso(const float value) {chargedIso_ = value;}
  void setneutralIso(const float value) {neutralIso_ = value;}
  void setphotonIso(const float value) {photonIso_ = value;}
  void setpuIso(const float value) {puIso_ = value;}
  void setdB(const float value) {dB_ = value;}
  void setipDXY(const float value) {ipDXY_ = value;}
  void setdz(const float value) {dz_ = value;}
  void setnMissingInnerHits(const float value) {nMissingInnerHits_ = value;}
  void setr9(const float value) {r9_ = value;}
  void setESCOverETrack(const float value) {ESCOverETrack_ = value;}
  void setDEtaSCTrk(const float value) {DEtaSCTrk_ = value;}
  void setDPhiSCTrk(const float value) {DPhiSCTrk_ = value;}
  void setecalEnergy(const float value) {ecalEnergy_ = value;}
  void setpassConversionVeto(const bool value) {passConversionVeto_ = value;}
  void setisEB(const bool value) {isEB_ = value;}
  void setisEE(const bool value) {isEE_ = value;}
  void setisEBGap(const bool value) {isEBGap_ = value;}
  void setisEBEtaGap(const bool value) {isEBEtaGap_ = value;}
  void setisEBPhiGap(const bool value) {isEBPhiGap_ = value;}
  void setisEEGap(const bool value) {isEEGap_ = value;}
  void setisEERingGap(const bool value) {isEERingGap_ = value;}
  void setisEEDeeGap(const bool value) {isEEDeeGap_ = value;}
  void setisEBEEGap(const bool value) {isEBEEGap_ = value;}
  void setisElectron(const bool value) {isElectron_ = value;}
  void setecalSeed(const bool value) {ecalSeed_ = value;}
  void settrackSeed(const bool value) {trackSeed_ = value;}
  void seteidCutLoose(const float value) {eidCutLoose_ = value;}
  void seteidCutMedium(const float value) {eidCutMedium_ = value;}
  void seteidCutTight(const float value) {eidCutTight_ = value;}
  void seteidCutVeto(const float value) {eidCutVeto_ = value;}
  void seteidMVAWP80(const float value) {eidMVAWP80_ = value;}
  void seteidMVAWP90(const float value) {eidMVAWP90_ = value;}
  void seteidTrgMVAWP80(const float value) {eidTrgMVAWP80_ = value;}
  void seteidTrgMVAWP90(const float value) {eidTrgMVAWP90_ = value;}
  void setpfHadronIso(const float value) {pfHadronIso_ = value;}
  void setpfNeutralIso(const float value) {pfNeutralIso_ = value;}
  void setpfPhotonIso(const float value) {pfPhotonIso_ = value;}
  void setHLT_Ele27_eta2p1_WPLoose_Gsf(const bool value) {HLT_Ele27_eta2p1_WPLoose_Gsf_ = value;}
  void setHLT_DoubleEle33_CaloIdL_GsfTrkIdVL(const bool value) {HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_ = value;}
  void setHLT_Ele22_eta2p1_WPLoose_Gsf(const bool value) {HLT_Ele22_eta2p1_WPLoose_Gsf_ = value;}
  void setHLT_Ele23_WPLoose_Gsf(const bool value) {HLT_Ele23_WPLoose_Gsf_ = value;}
  void setHLT_Ele27_WPLoose_Gsf(const bool value) {HLT_Ele27_WPLoose_Gsf_ = value;}
  void sete1x5(const float value) {e1x5_ = value;}
  void sete5x5(const float value) {e5x5_ = value;}
  void setsigmaIEtaIEta(const float value) {sigmaIEtaIEta_ = value;}
  void setfull5x5_sigmaIEtaIEta(const float value) {full5x5_sigmaIEtaIEta_ = value;}
  void setsigmaIPhiIPhi(const float value) {sigmaIPhiIPhi_ = value;}
  void sethadronicOverEM(const float value) {hadronicOverEM_ = value;}
  void setx(const float value) {x_ = value;}
  void sety(const float value) {y_ = value;}
  void setz(const float value) {z_ = value;}
  void setenergy(const float value) {energy_ = value;}
  void setrawEnergy(const float value) {rawEnergy_ = value;}
  void setphiWidth(const float value) {phiWidth_ = value;}
  void setetaWidth(const float value) {etaWidth_ = value;}
  void setLotentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Genjet: public TLorentzVector{
friend class URStreamer;
public:
//  Genjet(const int &i_charge_,const float &i_e_,const float &i_invisibleEnergy_,const float &i_pdgId_):
//    
//  {}
  Genjet():
    TLorentzVector(),
    charge_(0),
    e_(0),
    invisibleEnergy_(0),
    pdgId_(0)
  {}
  int charge() const {return charge_;}
  float e() const {return e_;}
  float invisibleEnergy() const {return invisibleEnergy_;}
  float pdgId() const {return pdgId_;}
private:
  int charge_;
  float e_;
  float invisibleEnergy_;
  float pdgId_;
  void setcharge(const int value) {charge_ = value;}
  void sete(const float value) {e_ = value;}
  void setinvisibleEnergy(const float value) {invisibleEnergy_ = value;}
  void setpdgId(const float value) {pdgId_ = value;}
  void setLotentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Rho{
friend class URStreamer;
public:
//  Rho(const Double_t &i_value_):
//    
//  {}
  Rho():
    value_(0)
  {}
  Double_t value() const {return value_;}
private:
  Double_t value_;
  void setvalue(const Double_t value) {value_ = value;}
};

class Jet: public TLorentzVector{
friend class URStreamer;
public:
//  Jet(const int &i_charge_,const float &i_e_,const float &i_area_,const float &i_JESUnc_,const float &i_JER_,const float &i_JERUp_,const float &i_JERDown_,const float &i_uncorrPt_,const float &i_uncorrEta_,const float &i_uncorrPhi_,const float &i_uncorrM_,const float &i_uncorrEnergy_,const float &i_chargedHadronEnergyFraction_,const float &i_neutralHadronEnergyFraction_,const float &i_chargedEmEnergyFraction_,const float &i_neutralEmEnergyFraction_,const float &i_HFHadronEnergyFraction_,const float &i_HFEMEnergyFraction_,const float &i_muonEnergyFraction_,const float &i_chargedMultiplicity_,const float &i_neutralMultiplicity_,const float &i_numChargedHadrons_,const float &i_numNeutralHadrons_,const float &i_numPhotons_,const float &i_numElectrons_,const float &i_numMuons_,const float &i_numForwardEMs_,const float &i_numForwardHads_,const float &i_numberOfDaughters_,const float &i_puId_,const float &i_jetBProb_,const float &i_jetProb_,const float &i_trkHiPur_,const float &i_trkHiEff_,const float &i_ssvHiEff_,const float &i_ssvHiPur_,const float &i_csv_,const float &i_csvIncl_,const float &i_CvsLtag_,const float &i_CombinedMVA_,const float &i_CvsBtag_,const float &i_vtxMass_,const float &i_vtxNtracks_,const float &i_vtx3DVal_,const float &i_vtx3DSig_,const int &i_partonFlavour_,const int &i_hadronFlavour_):
//    
//  {}
  Jet():
    TLorentzVector(),
    charge_(0),
    e_(0),
    area_(0),
    JESUnc_(0),
    JER_(0),
    JERUp_(0),
    JERDown_(0),
    uncorrPt_(0),
    uncorrEta_(0),
    uncorrPhi_(0),
    uncorrM_(0),
    uncorrEnergy_(0),
    chargedHadronEnergyFraction_(0),
    neutralHadronEnergyFraction_(0),
    chargedEmEnergyFraction_(0),
    neutralEmEnergyFraction_(0),
    HFHadronEnergyFraction_(0),
    HFEMEnergyFraction_(0),
    muonEnergyFraction_(0),
    chargedMultiplicity_(0),
    neutralMultiplicity_(0),
    numChargedHadrons_(0),
    numNeutralHadrons_(0),
    numPhotons_(0),
    numElectrons_(0),
    numMuons_(0),
    numForwardEMs_(0),
    numForwardHads_(0),
    numberOfDaughters_(0),
    puId_(0),
    jetBProb_(0),
    jetProb_(0),
    trkHiPur_(0),
    trkHiEff_(0),
    ssvHiEff_(0),
    ssvHiPur_(0),
    csv_(0),
    csvIncl_(0),
    CvsLtag_(0),
    CombinedMVA_(0),
    CvsBtag_(0),
    vtxMass_(0),
    vtxNtracks_(0),
    vtx3DVal_(0),
    vtx3DSig_(0),
    partonFlavour_(0),
    hadronFlavour_(0)
  {}
  int charge() const {return charge_;}
  float e() const {return e_;}
  float area() const {return area_;}
  float JESUnc() const {return JESUnc_;}
  float JER() const {return JER_;}
  float JERUp() const {return JERUp_;}
  float JERDown() const {return JERDown_;}
  float uncorrPt() const {return uncorrPt_;}
  float uncorrEta() const {return uncorrEta_;}
  float uncorrPhi() const {return uncorrPhi_;}
  float uncorrM() const {return uncorrM_;}
  float uncorrEnergy() const {return uncorrEnergy_;}
  float chargedHadronEnergyFraction() const {return chargedHadronEnergyFraction_;}
  float neutralHadronEnergyFraction() const {return neutralHadronEnergyFraction_;}
  float chargedEmEnergyFraction() const {return chargedEmEnergyFraction_;}
  float neutralEmEnergyFraction() const {return neutralEmEnergyFraction_;}
  float HFHadronEnergyFraction() const {return HFHadronEnergyFraction_;}
  float HFEMEnergyFraction() const {return HFEMEnergyFraction_;}
  float muonEnergyFraction() const {return muonEnergyFraction_;}
  float chargedMultiplicity() const {return chargedMultiplicity_;}
  float neutralMultiplicity() const {return neutralMultiplicity_;}
  float numChargedHadrons() const {return numChargedHadrons_;}
  float numNeutralHadrons() const {return numNeutralHadrons_;}
  float numPhotons() const {return numPhotons_;}
  float numElectrons() const {return numElectrons_;}
  float numMuons() const {return numMuons_;}
  float numForwardEMs() const {return numForwardEMs_;}
  float numForwardHads() const {return numForwardHads_;}
  float numberOfDaughters() const {return numberOfDaughters_;}
  float puId() const {return puId_;}
  float jetBProb() const {return jetBProb_;}
  float jetProb() const {return jetProb_;}
  float trkHiPur() const {return trkHiPur_;}
  float trkHiEff() const {return trkHiEff_;}
  float ssvHiEff() const {return ssvHiEff_;}
  float ssvHiPur() const {return ssvHiPur_;}
  float csv() const {return csv_;}
  float csvIncl() const {return csvIncl_;}
  float CvsLtag() const {return CvsLtag_;}
  float CombinedMVA() const {return CombinedMVA_;}
  float CvsBtag() const {return CvsBtag_;}
  float vtxMass() const {return vtxMass_;}
  float vtxNtracks() const {return vtxNtracks_;}
  float vtx3DVal() const {return vtx3DVal_;}
  float vtx3DSig() const {return vtx3DSig_;}
  int partonFlavour() const {return partonFlavour_;}
  int hadronFlavour() const {return hadronFlavour_;}
private:
  int charge_;
  float e_;
  float area_;
  float JESUnc_;
  float JER_;
  float JERUp_;
  float JERDown_;
  float uncorrPt_;
  float uncorrEta_;
  float uncorrPhi_;
  float uncorrM_;
  float uncorrEnergy_;
  float chargedHadronEnergyFraction_;
  float neutralHadronEnergyFraction_;
  float chargedEmEnergyFraction_;
  float neutralEmEnergyFraction_;
  float HFHadronEnergyFraction_;
  float HFEMEnergyFraction_;
  float muonEnergyFraction_;
  float chargedMultiplicity_;
  float neutralMultiplicity_;
  float numChargedHadrons_;
  float numNeutralHadrons_;
  float numPhotons_;
  float numElectrons_;
  float numMuons_;
  float numForwardEMs_;
  float numForwardHads_;
  float numberOfDaughters_;
  float puId_;
  float jetBProb_;
  float jetProb_;
  float trkHiPur_;
  float trkHiEff_;
  float ssvHiEff_;
  float ssvHiPur_;
  float csv_;
  float csvIncl_;
  float CvsLtag_;
  float CombinedMVA_;
  float CvsBtag_;
  float vtxMass_;
  float vtxNtracks_;
  float vtx3DVal_;
  float vtx3DSig_;
  int partonFlavour_;
  int hadronFlavour_;
  void setcharge(const int value) {charge_ = value;}
  void sete(const float value) {e_ = value;}
  void setarea(const float value) {area_ = value;}
  void setJESUnc(const float value) {JESUnc_ = value;}
  void setJER(const float value) {JER_ = value;}
  void setJERUp(const float value) {JERUp_ = value;}
  void setJERDown(const float value) {JERDown_ = value;}
  void setuncorrPt(const float value) {uncorrPt_ = value;}
  void setuncorrEta(const float value) {uncorrEta_ = value;}
  void setuncorrPhi(const float value) {uncorrPhi_ = value;}
  void setuncorrM(const float value) {uncorrM_ = value;}
  void setuncorrEnergy(const float value) {uncorrEnergy_ = value;}
  void setchargedHadronEnergyFraction(const float value) {chargedHadronEnergyFraction_ = value;}
  void setneutralHadronEnergyFraction(const float value) {neutralHadronEnergyFraction_ = value;}
  void setchargedEmEnergyFraction(const float value) {chargedEmEnergyFraction_ = value;}
  void setneutralEmEnergyFraction(const float value) {neutralEmEnergyFraction_ = value;}
  void setHFHadronEnergyFraction(const float value) {HFHadronEnergyFraction_ = value;}
  void setHFEMEnergyFraction(const float value) {HFEMEnergyFraction_ = value;}
  void setmuonEnergyFraction(const float value) {muonEnergyFraction_ = value;}
  void setchargedMultiplicity(const float value) {chargedMultiplicity_ = value;}
  void setneutralMultiplicity(const float value) {neutralMultiplicity_ = value;}
  void setnumChargedHadrons(const float value) {numChargedHadrons_ = value;}
  void setnumNeutralHadrons(const float value) {numNeutralHadrons_ = value;}
  void setnumPhotons(const float value) {numPhotons_ = value;}
  void setnumElectrons(const float value) {numElectrons_ = value;}
  void setnumMuons(const float value) {numMuons_ = value;}
  void setnumForwardEMs(const float value) {numForwardEMs_ = value;}
  void setnumForwardHads(const float value) {numForwardHads_ = value;}
  void setnumberOfDaughters(const float value) {numberOfDaughters_ = value;}
  void setpuId(const float value) {puId_ = value;}
  void setjetBProb(const float value) {jetBProb_ = value;}
  void setjetProb(const float value) {jetProb_ = value;}
  void settrkHiPur(const float value) {trkHiPur_ = value;}
  void settrkHiEff(const float value) {trkHiEff_ = value;}
  void setssvHiEff(const float value) {ssvHiEff_ = value;}
  void setssvHiPur(const float value) {ssvHiPur_ = value;}
  void setcsv(const float value) {csv_ = value;}
  void setcsvIncl(const float value) {csvIncl_ = value;}
  void setCvsLtag(const float value) {CvsLtag_ = value;}
  void setCombinedMVA(const float value) {CombinedMVA_ = value;}
  void setCvsBtag(const float value) {CvsBtag_ = value;}
  void setvtxMass(const float value) {vtxMass_ = value;}
  void setvtxNtracks(const float value) {vtxNtracks_ = value;}
  void setvtx3DVal(const float value) {vtx3DVal_ = value;}
  void setvtx3DSig(const float value) {vtx3DSig_ = value;}
  void setpartonFlavour(const int value) {partonFlavour_ = value;}
  void sethadronFlavour(const int value) {hadronFlavour_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Muon: public TLorentzVector{
friend class URStreamer;
public:
//  Muon(const int &i_charge_,const float &i_dB_,const float &i_ipDXY_,const float &i_dz_,const float &i_nMissingInnerHits_,const float &i_chargedIso_,const float &i_neutralIso_,const float &i_photonIso_,const float &i_puIso_,const float &i_ECalEnergy_,const float &i_HCalEnergy_,const int &i_numChambers_,const int &i_numMatchedStations_,const float &i_trackiso_,const float &i_ecaliso_,const float &i_hcaliso_,const float &i_pfChargedIso04_,const float &i_pfNeutralIso04_,const float &i_pfPhotonIso04_,const float &i_pfPUIso04_,const float &i_trkIso03_,const float &i_ptErr_,const float &i_chi2_,const int &i_ndof_,const float &i_validHits_,const float &i_pixelHits_,const float &i_trackerLayers_,const bool &i_isGlobal_,const bool &i_isTracker_,const bool &i_isCalo_,const bool &i_isPF_,const bool &i_isStandAlone_,const bool &i_isLoose_,const bool &i_HLT_DoubleIsoMu17_eta2p1_,const bool &i_HLT_IsoMu24_eta2p1_,const bool &i_HLT_IsoMu20_eta2p1_,const bool &i_HLT_IsoMu18_,const bool &i_HLT_IsoTkMu20_,const bool &i_HLT_IsoMu20_,const bool &i_HLT_IsoTkMu27_,const bool &i_HLT_IsoMu27_,const bool &i_HLT_Mu50_,const bool &i_HLT_Mu45_eta2p1_,const bool &i_HLT_IsoMu22_,const bool &i_HLT_IsoTkMu22_):
//    
//  {}
  Muon():
    TLorentzVector(),
    charge_(0),
    dB_(0),
    ipDXY_(0),
    dz_(0),
    nMissingInnerHits_(0),
    chargedIso_(0),
    neutralIso_(0),
    photonIso_(0),
    puIso_(0),
    ECalEnergy_(0),
    HCalEnergy_(0),
    numChambers_(0),
    numMatchedStations_(0),
    trackiso_(0),
    ecaliso_(0),
    hcaliso_(0),
    pfChargedIso04_(0),
    pfNeutralIso04_(0),
    pfPhotonIso04_(0),
    pfPUIso04_(0),
    trkIso03_(0),
    ptErr_(0),
    chi2_(0),
    ndof_(0),
    validHits_(0),
    pixelHits_(0),
    trackerLayers_(0),
    isGlobal_(0),
    isTracker_(0),
    isCalo_(0),
    isPF_(0),
    isStandAlone_(0),
    isLoose_(0),
    HLT_DoubleIsoMu17_eta2p1_(0),
    HLT_IsoMu24_eta2p1_(0),
    HLT_IsoMu20_eta2p1_(0),
    HLT_IsoMu18_(0),
    HLT_IsoTkMu20_(0),
    HLT_IsoMu20_(0),
    HLT_IsoTkMu27_(0),
    HLT_IsoMu27_(0),
    HLT_Mu50_(0),
    HLT_Mu45_eta2p1_(0),
    HLT_IsoMu22_(0),
    HLT_IsoTkMu22_(0)
  {}
  int charge() const {return charge_;}
  float dB() const {return dB_;}
  float ipDXY() const {return ipDXY_;}
  float dz() const {return dz_;}
  float nMissingInnerHits() const {return nMissingInnerHits_;}
  float chargedIso() const {return chargedIso_;}
  float neutralIso() const {return neutralIso_;}
  float photonIso() const {return photonIso_;}
  float puIso() const {return puIso_;}
  float ECalEnergy() const {return ECalEnergy_;}
  float HCalEnergy() const {return HCalEnergy_;}
  int numChambers() const {return numChambers_;}
  int numMatchedStations() const {return numMatchedStations_;}
  float trackiso() const {return trackiso_;}
  float ecaliso() const {return ecaliso_;}
  float hcaliso() const {return hcaliso_;}
  float pfChargedIso04() const {return pfChargedIso04_;}
  float pfNeutralIso04() const {return pfNeutralIso04_;}
  float pfPhotonIso04() const {return pfPhotonIso04_;}
  float pfPUIso04() const {return pfPUIso04_;}
  float trkIso03() const {return trkIso03_;}
  float ptErr() const {return ptErr_;}
  float chi2() const {return chi2_;}
  int ndof() const {return ndof_;}
  float validHits() const {return validHits_;}
  float pixelHits() const {return pixelHits_;}
  float trackerLayers() const {return trackerLayers_;}
  bool isGlobal() const {return isGlobal_;}
  bool isTracker() const {return isTracker_;}
  bool isCalo() const {return isCalo_;}
  bool isPF() const {return isPF_;}
  bool isStandAlone() const {return isStandAlone_;}
  bool isLoose() const {return isLoose_;}
  bool HLT_DoubleIsoMu17_eta2p1() const {return HLT_DoubleIsoMu17_eta2p1_;}
  bool HLT_IsoMu24_eta2p1() const {return HLT_IsoMu24_eta2p1_;}
  bool HLT_IsoMu20_eta2p1() const {return HLT_IsoMu20_eta2p1_;}
  bool HLT_IsoMu18() const {return HLT_IsoMu18_;}
  bool HLT_IsoTkMu20() const {return HLT_IsoTkMu20_;}
  bool HLT_IsoMu20() const {return HLT_IsoMu20_;}
  bool HLT_IsoTkMu27() const {return HLT_IsoTkMu27_;}
  bool HLT_IsoMu27() const {return HLT_IsoMu27_;}
  bool HLT_Mu50() const {return HLT_Mu50_;}
  bool HLT_Mu45_eta2p1() const {return HLT_Mu45_eta2p1_;}
  bool HLT_IsoMu22() const {return HLT_IsoMu22_;}
  bool HLT_IsoTkMu22() const {return HLT_IsoTkMu22_;}
private:
  int charge_;
  float dB_;
  float ipDXY_;
  float dz_;
  float nMissingInnerHits_;
  float chargedIso_;
  float neutralIso_;
  float photonIso_;
  float puIso_;
  float ECalEnergy_;
  float HCalEnergy_;
  int numChambers_;
  int numMatchedStations_;
  float trackiso_;
  float ecaliso_;
  float hcaliso_;
  float pfChargedIso04_;
  float pfNeutralIso04_;
  float pfPhotonIso04_;
  float pfPUIso04_;
  float trkIso03_;
  float ptErr_;
  float chi2_;
  int ndof_;
  float validHits_;
  float pixelHits_;
  float trackerLayers_;
  bool isGlobal_;
  bool isTracker_;
  bool isCalo_;
  bool isPF_;
  bool isStandAlone_;
  bool isLoose_;
  bool HLT_DoubleIsoMu17_eta2p1_;
  bool HLT_IsoMu24_eta2p1_;
  bool HLT_IsoMu20_eta2p1_;
  bool HLT_IsoMu18_;
  bool HLT_IsoTkMu20_;
  bool HLT_IsoMu20_;
  bool HLT_IsoTkMu27_;
  bool HLT_IsoMu27_;
  bool HLT_Mu50_;
  bool HLT_Mu45_eta2p1_;
  bool HLT_IsoMu22_;
  bool HLT_IsoTkMu22_;
  void setcharge(const int value) {charge_ = value;}
  void setdB(const float value) {dB_ = value;}
  void setipDXY(const float value) {ipDXY_ = value;}
  void setdz(const float value) {dz_ = value;}
  void setnMissingInnerHits(const float value) {nMissingInnerHits_ = value;}
  void setchargedIso(const float value) {chargedIso_ = value;}
  void setneutralIso(const float value) {neutralIso_ = value;}
  void setphotonIso(const float value) {photonIso_ = value;}
  void setpuIso(const float value) {puIso_ = value;}
  void setECalEnergy(const float value) {ECalEnergy_ = value;}
  void setHCalEnergy(const float value) {HCalEnergy_ = value;}
  void setnumChambers(const int value) {numChambers_ = value;}
  void setnumMatchedStations(const int value) {numMatchedStations_ = value;}
  void settrackiso(const float value) {trackiso_ = value;}
  void setecaliso(const float value) {ecaliso_ = value;}
  void sethcaliso(const float value) {hcaliso_ = value;}
  void setpfChargedIso04(const float value) {pfChargedIso04_ = value;}
  void setpfNeutralIso04(const float value) {pfNeutralIso04_ = value;}
  void setpfPhotonIso04(const float value) {pfPhotonIso04_ = value;}
  void setpfPUIso04(const float value) {pfPUIso04_ = value;}
  void settrkIso03(const float value) {trkIso03_ = value;}
  void setptErr(const float value) {ptErr_ = value;}
  void setchi2(const float value) {chi2_ = value;}
  void setndof(const int value) {ndof_ = value;}
  void setvalidHits(const float value) {validHits_ = value;}
  void setpixelHits(const float value) {pixelHits_ = value;}
  void settrackerLayers(const float value) {trackerLayers_ = value;}
  void setisGlobal(const bool value) {isGlobal_ = value;}
  void setisTracker(const bool value) {isTracker_ = value;}
  void setisCalo(const bool value) {isCalo_ = value;}
  void setisPF(const bool value) {isPF_ = value;}
  void setisStandAlone(const bool value) {isStandAlone_ = value;}
  void setisLoose(const bool value) {isLoose_ = value;}
  void setHLT_DoubleIsoMu17_eta2p1(const bool value) {HLT_DoubleIsoMu17_eta2p1_ = value;}
  void setHLT_IsoMu24_eta2p1(const bool value) {HLT_IsoMu24_eta2p1_ = value;}
  void setHLT_IsoMu20_eta2p1(const bool value) {HLT_IsoMu20_eta2p1_ = value;}
  void setHLT_IsoMu18(const bool value) {HLT_IsoMu18_ = value;}
  void setHLT_IsoTkMu20(const bool value) {HLT_IsoTkMu20_ = value;}
  void setHLT_IsoMu20(const bool value) {HLT_IsoMu20_ = value;}
  void setHLT_IsoTkMu27(const bool value) {HLT_IsoTkMu27_ = value;}
  void setHLT_IsoMu27(const bool value) {HLT_IsoMu27_ = value;}
  void setHLT_Mu50(const bool value) {HLT_Mu50_ = value;}
  void setHLT_Mu45_eta2p1(const bool value) {HLT_Mu45_eta2p1_ = value;}
  void setHLT_IsoMu22(const bool value) {HLT_IsoMu22_ = value;}
  void setHLT_IsoTkMu22(const bool value) {HLT_IsoTkMu22_ = value;}
  void setLotentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Genparticle: public TLorentzVector{
friend class URStreamer;
public:
//  Genparticle(const int &i_charge_,const float &i_e_,const float &i_vx_,const float &i_vy_,const float &i_vz_,const int &i_pdgId_,const int &i_status_,const int &i_idx_,const vector<int> &i_momIdx_,const int &i_nDaught_,const int &i_firstDaughtIdx_):
//    
//  {}
  Genparticle():
    TLorentzVector(),
    charge_(0),
    e_(0),
    vx_(0),
    vy_(0),
    vz_(0),
    pdgId_(0),
    status_(0),
    idx_(0),
    momIdx_(0),
    nDaught_(0),
    firstDaughtIdx_(0)
  {}
  int charge() const {return charge_;}
  float e() const {return e_;}
  float vx() const {return vx_;}
  float vy() const {return vy_;}
  float vz() const {return vz_;}
  int pdgId() const {return pdgId_;}
  int status() const {return status_;}
  int idx() const {return idx_;}
  vector<int> momIdx() const {return momIdx_;}
  int nDaught() const {return nDaught_;}
  int firstDaughtIdx() const {return firstDaughtIdx_;}
private:
  int charge_;
  float e_;
  float vx_;
  float vy_;
  float vz_;
  int pdgId_;
  int status_;
  int idx_;
  vector<int> momIdx_;
  int nDaught_;
  int firstDaughtIdx_;
  void setcharge(const int value) {charge_ = value;}
  void sete(const float value) {e_ = value;}
  void setvx(const float value) {vx_ = value;}
  void setvy(const float value) {vy_ = value;}
  void setvz(const float value) {vz_ = value;}
  void setpdgId(const int value) {pdgId_ = value;}
  void setstatus(const int value) {status_ = value;}
  void setidx(const int value) {idx_ = value;}
  void setmomIdx(const vector<int> value) {momIdx_ = value;}
  void setnDaught(const int value) {nDaught_ = value;}
  void setfirstDaughtIdx(const int value) {firstDaughtIdx_ = value;}
  void setLotentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Vertex{
friend class URStreamer;
public:
//  Vertex(const float &i_x_,const float &i_y_,const float &i_z_,const float &i_chi2_,const float &i_ndof_,const float &i_nTracks_):
//    
//  {}
  Vertex():
    x_(0),
    y_(0),
    z_(0),
    chi2_(0),
    ndof_(0),
    nTracks_(0)
  {}
  float x() const {return x_;}
  float y() const {return y_;}
  float z() const {return z_;}
  float chi2() const {return chi2_;}
  float ndof() const {return ndof_;}
  float nTracks() const {return nTracks_;}
private:
  float x_;
  float y_;
  float z_;
  float chi2_;
  float ndof_;
  float nTracks_;
  void setx(const float value) {x_ = value;}
  void sety(const float value) {y_ = value;}
  void setz(const float value) {z_ = value;}
  void setchi2(const float value) {chi2_ = value;}
  void setndof(const float value) {ndof_ = value;}
  void setnTracks(const float value) {nTracks_ = value;}
};

class Puinfo{
friend class URStreamer;
public:
//  Puinfo(const float &i_bx_,const float &i_nPU_,const float &i_nInteractions_):
//    
//  {}
  Puinfo():
    bx_(0),
    nPU_(0),
    nInteractions_(0)
  {}
  float bx() const {return bx_;}
  float nPU() const {return nPU_;}
  float nInteractions() const {return nInteractions_;}
private:
  float bx_;
  float nPU_;
  float nInteractions_;
  void setbx(const float value) {bx_ = value;}
  void setnPU(const float value) {nPU_ = value;}
  void setnInteractions(const float value) {nInteractions_ = value;}
};


class URStreamer{
public:
  UInt_t run;
  UInt_t lumi;
  UInt_t evt;

  URStreamer(TTree *tree):
    run(0),
    lumi(0),
    evt(0),
    trigger_HLT_Ele27_eta2p1_WPLoose_Gsf_(0),
    trigger_HLT_DoubleIsoMu17_eta2p1_(0),
    trigger_HLT_notexists_(0),
    trigger_HLT_IsoMu24_eta2p1_(0),
    trigger_HLT_IsoMu20_eta2p1_(0),
    trigger_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_(0),
    trigger_HLT_IsoMu18_(0),
    trigger_HLT_IsoTkMu20_(0),
    trigger_HLT_IsoMu20_(0),
    trigger_HLT_IsoTkMu27_(0),
    trigger_HLT_IsoMu27_(0),
    trigger_HLT_Mu50_(0),
    trigger_HLT_Mu45_eta2p1_(0),
    trigger_HLT_Ele22_eta2p1_WPLoose_Gsf_(0),
    trigger_HLT_Ele23_WPLoose_Gsf_(0),
    trigger_HLT_Ele27_WPLoose_Gsf_(0),
    trigger_HLT_IsoMu22_(0),
    trigger_HLT_IsoTkMu22_(0),
    filter_Flag_goodVertices_(0),
    filter_Flag_CSCTightHaloFilter_(0),
    filter_Flag_trkPOGFilters_(0),
    filter_Flag_trkPOG_logErrorTooManyClusters_(0),
    filter_Flag_EcalDeadCellTriggerPrimitiveFilter_(0),
    filter_Flag_ecalLaserCorrFilter_(0),
    filter_Flag_trkPOG_manystripclus53X_(0),
    filter_Flag_eeBadScFilter_(0),
    filter_Flag_METFilters_(0),
    filter_Flag_HBHENoiseFilter_(0),
    filter_Flag_trkPOG_toomanystripclus53X_(0),
    filter_Flag_hcalLaserEventFilter_(0),
    filter_Flag_HBHENoiseIsoFilter_(0),
    filter_Flag_CSCTightHalo2015Filter_(0),
    rho_value_(0),
    muons_pt_(0),
    muons_eta_(0),
    muons_phi_(0),
    muons_charge_(0),
    muons_dB_(0),
    muons_ipDXY_(0),
    muons_dz_(0),
    muons_nMissingInnerHits_(0),
    muons_chargedIso_(0),
    muons_neutralIso_(0),
    muons_photonIso_(0),
    muons_puIso_(0),
    muons_ECalEnergy_(0),
    muons_HCalEnergy_(0),
    muons_numChambers_(0),
    muons_numMatchedStations_(0),
    muons_trackiso_(0),
    muons_ecaliso_(0),
    muons_hcaliso_(0),
    muons_pfChargedIso04_(0),
    muons_pfNeutralIso04_(0),
    muons_pfPhotonIso04_(0),
    muons_pfPUIso04_(0),
    muons_trkIso03_(0),
    muons_ptErr_(0),
    muons_chi2_(0),
    muons_ndof_(0),
    muons_validHits_(0),
    muons_pixelHits_(0),
    muons_trackerLayers_(0),
    muons_isGlobal_(0),
    muons_isTracker_(0),
    muons_isCalo_(0),
    muons_isPF_(0),
    muons_isStandAlone_(0),
    muons_isLoose_(0),
    muons_HLT_DoubleIsoMu17_eta2p1_(0),
    muons_HLT_IsoMu24_eta2p1_(0),
    muons_HLT_IsoMu20_eta2p1_(0),
    muons_HLT_IsoMu18_(0),
    muons_HLT_IsoTkMu20_(0),
    muons_HLT_IsoMu20_(0),
    muons_HLT_IsoTkMu27_(0),
    muons_HLT_IsoMu27_(0),
    muons_HLT_Mu50_(0),
    muons_HLT_Mu45_eta2p1_(0),
    muons_HLT_IsoMu22_(0),
    muons_HLT_IsoTkMu22_(0),
    jets_pt_(0),
    jets_eta_(0),
    jets_phi_(0),
    jets_charge_(0),
    jets_e_(0),
    jets_area_(0),
    jets_mass_(0),
    jets_JESUnc_(0),
    jets_JER_(0),
    jets_JERUp_(0),
    jets_JERDown_(0),
    jets_uncorrPt_(0),
    jets_uncorrEta_(0),
    jets_uncorrPhi_(0),
    jets_uncorrM_(0),
    jets_uncorrEnergy_(0),
    jets_chargedHadronEnergyFraction_(0),
    jets_neutralHadronEnergyFraction_(0),
    jets_chargedEmEnergyFraction_(0),
    jets_neutralEmEnergyFraction_(0),
    jets_HFHadronEnergyFraction_(0),
    jets_HFEMEnergyFraction_(0),
    jets_muonEnergyFraction_(0),
    jets_chargedMultiplicity_(0),
    jets_neutralMultiplicity_(0),
    jets_numChargedHadrons_(0),
    jets_numNeutralHadrons_(0),
    jets_numPhotons_(0),
    jets_numElectrons_(0),
    jets_numMuons_(0),
    jets_numForwardEMs_(0),
    jets_numForwardHads_(0),
    jets_numberOfDaughters_(0),
    jets_puId_(0),
    jets_jetBProb_(0),
    jets_jetProb_(0),
    jets_trkHiPur_(0),
    jets_trkHiEff_(0),
    jets_ssvHiEff_(0),
    jets_ssvHiPur_(0),
    jets_csv_(0),
    jets_csvIncl_(0),
    jets_CvsLtag_(0),
    jets_CombinedMVA_(0),
    jets_CvsBtag_(0),
    jets_vtxMass_(0),
    jets_vtxNtracks_(0),
    jets_vtx3DVal_(0),
    jets_vtx3DSig_(0),
    jets_partonFlavour_(0),
    jets_hadronFlavour_(0),
    electrons_pt_(0),
    electrons_eta_(0),
    electrons_phi_(0),
    electrons_charge_(0),
    electrons_chargedIso_(0),
    electrons_neutralIso_(0),
    electrons_photonIso_(0),
    electrons_puIso_(0),
    electrons_dB_(0),
    electrons_ipDXY_(0),
    electrons_dz_(0),
    electrons_nMissingInnerHits_(0),
    electrons_r9_(0),
    electrons_ESCOverETrack_(0),
    electrons_DEtaSCTrk_(0),
    electrons_DPhiSCTrk_(0),
    electrons_ecalEnergy_(0),
    electrons_passConversionVeto_(0),
    electrons_isEB_(0),
    electrons_isEE_(0),
    electrons_isEBGap_(0),
    electrons_isEBEtaGap_(0),
    electrons_isEBPhiGap_(0),
    electrons_isEEGap_(0),
    electrons_isEERingGap_(0),
    electrons_isEEDeeGap_(0),
    electrons_isEBEEGap_(0),
    electrons_isElectron_(0),
    electrons_ecalSeed_(0),
    electrons_trackSeed_(0),
    electrons_eidCutLoose_(0),
    electrons_eidCutMedium_(0),
    electrons_eidCutTight_(0),
    electrons_eidCutVeto_(0),
    electrons_eidMVAWP80_(0),
    electrons_eidMVAWP90_(0),
    electrons_eidTrgMVAWP80_(0),
    electrons_eidTrgMVAWP90_(0),
    electrons_pfHadronIso_(0),
    electrons_pfNeutralIso_(0),
    electrons_pfPhotonIso_(0),
    electrons_HLT_Ele27_eta2p1_WPLoose_Gsf_(0),
    electrons_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_(0),
    electrons_HLT_Ele22_eta2p1_WPLoose_Gsf_(0),
    electrons_HLT_Ele23_WPLoose_Gsf_(0),
    electrons_HLT_Ele27_WPLoose_Gsf_(0),
    electrons_e1x5_(0),
    electrons_e5x5_(0),
    electrons_sigmaIEtaIEta_(0),
    electrons_full5x5_sigmaIEtaIEta_(0),
    electrons_sigmaIPhiIPhi_(0),
    electrons_hadronicOverEM_(0),
    electrons_x_(0),
    electrons_y_(0),
    electrons_z_(0),
    electrons_energy_(0),
    electrons_rawEnergy_(0),
    electrons_phiWidth_(0),
    electrons_etaWidth_(0),
    photons_pt_(0),
    photons_eta_(0),
    photons_phi_(0),
    photons_charge_(0),
    photons_x_(0),
    photons_y_(0),
    photons_z_(0),
    photons_energy_(0),
    photons_rawEnergy_(0),
    photons_phiWidth_(0),
    photons_etaWidth_(0),
    photons_e3x3_(0),
    photons_maxCrystalEnergy_(0),
    photons_isEB_(0),
    photons_isEE_(0),
    photons_isPFlowPhoton_(0),
    photons_hasConversionTracks_(0),
    photons_hasPixelSeed_(0),
    vertexs_x_(0),
    vertexs_y_(0),
    vertexs_z_(0),
    vertexs_chi2_(0),
    vertexs_ndof_(0),
    vertexs_nTracks_(0),
    METs_px_(0),
    METs_py_(0),
    METs_pxsmear_(0),
    METs_pysmear_(0),
    METs_pxunc_(0),
    METs_pyunc_(0),
    METs_pxuncJES_(0),
    METs_pyuncJES_(0),
    METs_pxuncJER_(0),
    METs_pyuncJER_(0),
    genInfo_weight_(0),
    genInfo_pdfid1_(0),
    genInfo_pdfid2_(0),
    genInfo_x1_(0),
    genInfo_x2_(0),
    genInfo_renScale_(0),
    MCWeights_weights_(0),
    NPNLOLHE_npnlo_(0),
    PUInfos_bx_(0),
    PUInfos_nPU_(0),
    PUInfos_nInteractions_(0),
    genParticles_pt_(0),
    genParticles_eta_(0),
    genParticles_phi_(0),
    genParticles_charge_(0),
    genParticles_e_(0),
    genParticles_vx_(0),
    genParticles_vy_(0),
    genParticles_vz_(0),
    genParticles_pdgId_(0),
    genParticles_status_(0),
    LHEPaticles_px_(0),
    LHEPaticles_py_(0),
    LHEPaticles_pz_(0),
    LHEPaticles_e_(0),
    LHEPaticles_pdgid_(0),
    LHEPaticles_status_(0),
    LHEPaticles_fmother_(0),
    LHEPaticles_lmother_(0),
    genParticles_idx_(0),
    genParticles_momIdx_(0),
    genParticles_nDaught_(0),
    genParticles_firstDaughtIdx_(0),
    genjets_pt_(0),
    genjets_eta_(0),
    genjets_phi_(0),
    genjets_charge_(0),
    genjets_e_(0),
    genjets_invisibleEnergy_(0),
    genjets_pdgId_(0),
    are_NPNLOLHE_loaded_(0), NPNLOLHE_(),
    are_METs_loaded_(0), METs_(),
    are_MCWeights_loaded_(0), MCWeights_(),
    are_genInfo_loaded_(0), genInfo_(),
    are_photons_loaded_(0), photons_(),
    are_LHEPaticles_loaded_(0), LHEPaticles_(),
    are_filter_loaded_(0), filter_(),
    are_trigger_loaded_(0), trigger_(),
    are_electrons_loaded_(0), electrons_(),
    are_genjets_loaded_(0), genjets_(),
    are_rho_loaded_(0), rho_(),
    are_jets_loaded_(0), jets_(),
    are_muons_loaded_(0), muons_(),
    are_genParticles_loaded_(0), genParticles_(),
    are_vertexs_loaded_(0), vertexs_(),
    are_PUInfos_loaded_(0), PUInfos_()
  {
    tree_ = tree;
    current_entry_ = -1;
    entries_ = tree_->GetEntries();
    tree_->SetBranchStatus("*",0); 
    tree_->SetBranchStatus("run", 1); tree_->SetBranchAddress("run", &run);
    tree_->SetBranchStatus("lumi", 1); tree_->SetBranchAddress("lumi", &lumi);
    tree_->SetBranchStatus("evt", 1); tree_->SetBranchAddress("evt", &evt);
  }

  ~URStreamer()
  {
    //{ EVT_DESTROY }
  }

  bool next(){
    
    METs_.clear();
    MCWeights_.clear();
    
    photons_.clear();
    LHEPaticles_.clear();
    
    
    electrons_.clear();
    genjets_.clear();
    
    jets_.clear();
    muons_.clear();
    genParticles_.clear();
    vertexs_.clear();
    PUInfos_.clear();
    current_entry_++;
    if(current_entry_ < entries_){
      tree_->GetEntry(current_entry_);
      return true;
    }
    return false;
  }

  void loadNpnlolhe(){
    if(!are_NPNLOLHE_loaded_){
      tree_->SetBranchStatus("NPNLOLHE.npnlo", 1); tree_->SetBranchAddress("NPNLOLHE.npnlo", &NPNLOLHE_npnlo_);
      are_NPNLOLHE_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadMets(){
    if(!are_METs_loaded_){
      tree_->SetBranchStatus("METs.px", 1); tree_->SetBranchAddress("METs.px", &METs_px_);
      tree_->SetBranchStatus("METs.py", 1); tree_->SetBranchAddress("METs.py", &METs_py_);
      tree_->SetBranchStatus("METs.pxsmear", 1); tree_->SetBranchAddress("METs.pxsmear", &METs_pxsmear_);
      tree_->SetBranchStatus("METs.pysmear", 1); tree_->SetBranchAddress("METs.pysmear", &METs_pysmear_);
      tree_->SetBranchStatus("METs.pxunc", 1); tree_->SetBranchAddress("METs.pxunc", &METs_pxunc_);
      tree_->SetBranchStatus("METs.pyunc", 1); tree_->SetBranchAddress("METs.pyunc", &METs_pyunc_);
      tree_->SetBranchStatus("METs.pxuncJES", 1); tree_->SetBranchAddress("METs.pxuncJES", &METs_pxuncJES_);
      tree_->SetBranchStatus("METs.pyuncJES", 1); tree_->SetBranchAddress("METs.pyuncJES", &METs_pyuncJES_);
      tree_->SetBranchStatus("METs.pxuncJER", 1); tree_->SetBranchAddress("METs.pxuncJER", &METs_pxuncJER_);
      tree_->SetBranchStatus("METs.pyuncJER", 1); tree_->SetBranchAddress("METs.pyuncJER", &METs_pyuncJER_);
      are_METs_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadMcweights(){
    if(!are_MCWeights_loaded_){
      tree_->SetBranchStatus("MCWeights.weights", 1); tree_->SetBranchAddress("MCWeights.weights", &MCWeights_weights_);
      are_MCWeights_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGeninfo(){
    if(!are_genInfo_loaded_){
      tree_->SetBranchStatus("genInfo.weight", 1); tree_->SetBranchAddress("genInfo.weight", &genInfo_weight_);
      tree_->SetBranchStatus("genInfo.pdfid1", 1); tree_->SetBranchAddress("genInfo.pdfid1", &genInfo_pdfid1_);
      tree_->SetBranchStatus("genInfo.pdfid2", 1); tree_->SetBranchAddress("genInfo.pdfid2", &genInfo_pdfid2_);
      tree_->SetBranchStatus("genInfo.x1", 1); tree_->SetBranchAddress("genInfo.x1", &genInfo_x1_);
      tree_->SetBranchStatus("genInfo.x2", 1); tree_->SetBranchAddress("genInfo.x2", &genInfo_x2_);
      tree_->SetBranchStatus("genInfo.renScale", 1); tree_->SetBranchAddress("genInfo.renScale", &genInfo_renScale_);
      are_genInfo_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadPhotons(){
    if(!are_photons_loaded_){
      tree_->SetBranchStatus("photons.pt", 1); tree_->SetBranchAddress("photons.pt", &photons_pt_);
      tree_->SetBranchStatus("photons.eta", 1); tree_->SetBranchAddress("photons.eta", &photons_eta_);
      tree_->SetBranchStatus("photons.phi", 1); tree_->SetBranchAddress("photons.phi", &photons_phi_);
      tree_->SetBranchStatus("photons.charge", 1); tree_->SetBranchAddress("photons.charge", &photons_charge_);
      tree_->SetBranchStatus("photons.x", 1); tree_->SetBranchAddress("photons.x", &photons_x_);
      tree_->SetBranchStatus("photons.y", 1); tree_->SetBranchAddress("photons.y", &photons_y_);
      tree_->SetBranchStatus("photons.z", 1); tree_->SetBranchAddress("photons.z", &photons_z_);
      tree_->SetBranchStatus("photons.energy", 1); tree_->SetBranchAddress("photons.energy", &photons_energy_);
      tree_->SetBranchStatus("photons.rawEnergy", 1); tree_->SetBranchAddress("photons.rawEnergy", &photons_rawEnergy_);
      tree_->SetBranchStatus("photons.phiWidth", 1); tree_->SetBranchAddress("photons.phiWidth", &photons_phiWidth_);
      tree_->SetBranchStatus("photons.etaWidth", 1); tree_->SetBranchAddress("photons.etaWidth", &photons_etaWidth_);
      tree_->SetBranchStatus("photons.e3x3", 1); tree_->SetBranchAddress("photons.e3x3", &photons_e3x3_);
      tree_->SetBranchStatus("photons.maxCrystalEnergy", 1); tree_->SetBranchAddress("photons.maxCrystalEnergy", &photons_maxCrystalEnergy_);
      tree_->SetBranchStatus("photons.isEB", 1); tree_->SetBranchAddress("photons.isEB", &photons_isEB_);
      tree_->SetBranchStatus("photons.isEE", 1); tree_->SetBranchAddress("photons.isEE", &photons_isEE_);
      tree_->SetBranchStatus("photons.isPFlowPhoton", 1); tree_->SetBranchAddress("photons.isPFlowPhoton", &photons_isPFlowPhoton_);
      tree_->SetBranchStatus("photons.hasConversionTracks", 1); tree_->SetBranchAddress("photons.hasConversionTracks", &photons_hasConversionTracks_);
      tree_->SetBranchStatus("photons.hasPixelSeed", 1); tree_->SetBranchAddress("photons.hasPixelSeed", &photons_hasPixelSeed_);
      are_photons_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadLhepaticles(){
    if(!are_LHEPaticles_loaded_){
      tree_->SetBranchStatus("LHEPaticles.px", 1); tree_->SetBranchAddress("LHEPaticles.px", &LHEPaticles_px_);
      tree_->SetBranchStatus("LHEPaticles.py", 1); tree_->SetBranchAddress("LHEPaticles.py", &LHEPaticles_py_);
      tree_->SetBranchStatus("LHEPaticles.pz", 1); tree_->SetBranchAddress("LHEPaticles.pz", &LHEPaticles_pz_);
      tree_->SetBranchStatus("LHEPaticles.e", 1); tree_->SetBranchAddress("LHEPaticles.e", &LHEPaticles_e_);
      tree_->SetBranchStatus("LHEPaticles.pdgid", 1); tree_->SetBranchAddress("LHEPaticles.pdgid", &LHEPaticles_pdgid_);
      tree_->SetBranchStatus("LHEPaticles.status", 1); tree_->SetBranchAddress("LHEPaticles.status", &LHEPaticles_status_);
      tree_->SetBranchStatus("LHEPaticles.fmother", 1); tree_->SetBranchAddress("LHEPaticles.fmother", &LHEPaticles_fmother_);
      tree_->SetBranchStatus("LHEPaticles.lmother", 1); tree_->SetBranchAddress("LHEPaticles.lmother", &LHEPaticles_lmother_);
      are_LHEPaticles_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadFilter(){
    if(!are_filter_loaded_){
      tree_->SetBranchStatus("filter.Flag_goodVertices", 1); tree_->SetBranchAddress("filter.Flag_goodVertices", &filter_Flag_goodVertices_);
      tree_->SetBranchStatus("filter.Flag_CSCTightHaloFilter", 1); tree_->SetBranchAddress("filter.Flag_CSCTightHaloFilter", &filter_Flag_CSCTightHaloFilter_);
      tree_->SetBranchStatus("filter.Flag_trkPOGFilters", 1); tree_->SetBranchAddress("filter.Flag_trkPOGFilters", &filter_Flag_trkPOGFilters_);
      tree_->SetBranchStatus("filter.Flag_trkPOG_logErrorTooManyClusters", 1); tree_->SetBranchAddress("filter.Flag_trkPOG_logErrorTooManyClusters", &filter_Flag_trkPOG_logErrorTooManyClusters_);
      tree_->SetBranchStatus("filter.Flag_EcalDeadCellTriggerPrimitiveFilter", 1); tree_->SetBranchAddress("filter.Flag_EcalDeadCellTriggerPrimitiveFilter", &filter_Flag_EcalDeadCellTriggerPrimitiveFilter_);
      tree_->SetBranchStatus("filter.Flag_ecalLaserCorrFilter", 1); tree_->SetBranchAddress("filter.Flag_ecalLaserCorrFilter", &filter_Flag_ecalLaserCorrFilter_);
      tree_->SetBranchStatus("filter.Flag_trkPOG_manystripclus53X", 1); tree_->SetBranchAddress("filter.Flag_trkPOG_manystripclus53X", &filter_Flag_trkPOG_manystripclus53X_);
      tree_->SetBranchStatus("filter.Flag_eeBadScFilter", 1); tree_->SetBranchAddress("filter.Flag_eeBadScFilter", &filter_Flag_eeBadScFilter_);
      tree_->SetBranchStatus("filter.Flag_METFilters", 1); tree_->SetBranchAddress("filter.Flag_METFilters", &filter_Flag_METFilters_);
      tree_->SetBranchStatus("filter.Flag_HBHENoiseFilter", 1); tree_->SetBranchAddress("filter.Flag_HBHENoiseFilter", &filter_Flag_HBHENoiseFilter_);
      tree_->SetBranchStatus("filter.Flag_trkPOG_toomanystripclus53X", 1); tree_->SetBranchAddress("filter.Flag_trkPOG_toomanystripclus53X", &filter_Flag_trkPOG_toomanystripclus53X_);
      tree_->SetBranchStatus("filter.Flag_hcalLaserEventFilter", 1); tree_->SetBranchAddress("filter.Flag_hcalLaserEventFilter", &filter_Flag_hcalLaserEventFilter_);
      tree_->SetBranchStatus("filter.Flag_HBHENoiseIsoFilter", 1); tree_->SetBranchAddress("filter.Flag_HBHENoiseIsoFilter", &filter_Flag_HBHENoiseIsoFilter_);
      tree_->SetBranchStatus("filter.Flag_CSCTightHalo2015Filter", 1); tree_->SetBranchAddress("filter.Flag_CSCTightHalo2015Filter", &filter_Flag_CSCTightHalo2015Filter_);
      are_filter_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadTrigger(){
    if(!are_trigger_loaded_){
      tree_->SetBranchStatus("trigger.HLT_Ele27_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("trigger.HLT_Ele27_eta2p1_WPLoose_Gsf", &trigger_HLT_Ele27_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("trigger.HLT_DoubleIsoMu17_eta2p1", 1); tree_->SetBranchAddress("trigger.HLT_DoubleIsoMu17_eta2p1", &trigger_HLT_DoubleIsoMu17_eta2p1_);
      tree_->SetBranchStatus("trigger.HLT_notexists", 1); tree_->SetBranchAddress("trigger.HLT_notexists", &trigger_HLT_notexists_);
      tree_->SetBranchStatus("trigger.HLT_IsoMu24_eta2p1", 1); tree_->SetBranchAddress("trigger.HLT_IsoMu24_eta2p1", &trigger_HLT_IsoMu24_eta2p1_);
      tree_->SetBranchStatus("trigger.HLT_IsoMu20_eta2p1", 1); tree_->SetBranchAddress("trigger.HLT_IsoMu20_eta2p1", &trigger_HLT_IsoMu20_eta2p1_);
      tree_->SetBranchStatus("trigger.HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("trigger.HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &trigger_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("trigger.HLT_IsoMu18", 1); tree_->SetBranchAddress("trigger.HLT_IsoMu18", &trigger_HLT_IsoMu18_);
      tree_->SetBranchStatus("trigger.HLT_IsoTkMu20", 1); tree_->SetBranchAddress("trigger.HLT_IsoTkMu20", &trigger_HLT_IsoTkMu20_);
      tree_->SetBranchStatus("trigger.HLT_IsoMu20", 1); tree_->SetBranchAddress("trigger.HLT_IsoMu20", &trigger_HLT_IsoMu20_);
      tree_->SetBranchStatus("trigger.HLT_IsoTkMu27", 1); tree_->SetBranchAddress("trigger.HLT_IsoTkMu27", &trigger_HLT_IsoTkMu27_);
      tree_->SetBranchStatus("trigger.HLT_IsoMu27", 1); tree_->SetBranchAddress("trigger.HLT_IsoMu27", &trigger_HLT_IsoMu27_);
      tree_->SetBranchStatus("trigger.HLT_Mu50", 1); tree_->SetBranchAddress("trigger.HLT_Mu50", &trigger_HLT_Mu50_);
      tree_->SetBranchStatus("trigger.HLT_Mu45_eta2p1", 1); tree_->SetBranchAddress("trigger.HLT_Mu45_eta2p1", &trigger_HLT_Mu45_eta2p1_);
      tree_->SetBranchStatus("trigger.HLT_Ele22_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("trigger.HLT_Ele22_eta2p1_WPLoose_Gsf", &trigger_HLT_Ele22_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("trigger.HLT_Ele23_WPLoose_Gsf", 1); tree_->SetBranchAddress("trigger.HLT_Ele23_WPLoose_Gsf", &trigger_HLT_Ele23_WPLoose_Gsf_);
      tree_->SetBranchStatus("trigger.HLT_Ele27_WPLoose_Gsf", 1); tree_->SetBranchAddress("trigger.HLT_Ele27_WPLoose_Gsf", &trigger_HLT_Ele27_WPLoose_Gsf_);
      tree_->SetBranchStatus("trigger.HLT_IsoMu22", 1); tree_->SetBranchAddress("trigger.HLT_IsoMu22", &trigger_HLT_IsoMu22_);
      tree_->SetBranchStatus("trigger.HLT_IsoTkMu22", 1); tree_->SetBranchAddress("trigger.HLT_IsoTkMu22", &trigger_HLT_IsoTkMu22_);
      are_trigger_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadElectrons(){
    if(!are_electrons_loaded_){
      tree_->SetBranchStatus("electrons.pt", 1); tree_->SetBranchAddress("electrons.pt", &electrons_pt_);
      tree_->SetBranchStatus("electrons.eta", 1); tree_->SetBranchAddress("electrons.eta", &electrons_eta_);
      tree_->SetBranchStatus("electrons.phi", 1); tree_->SetBranchAddress("electrons.phi", &electrons_phi_);
      tree_->SetBranchStatus("electrons.charge", 1); tree_->SetBranchAddress("electrons.charge", &electrons_charge_);
      tree_->SetBranchStatus("electrons.chargedIso", 1); tree_->SetBranchAddress("electrons.chargedIso", &electrons_chargedIso_);
      tree_->SetBranchStatus("electrons.neutralIso", 1); tree_->SetBranchAddress("electrons.neutralIso", &electrons_neutralIso_);
      tree_->SetBranchStatus("electrons.photonIso", 1); tree_->SetBranchAddress("electrons.photonIso", &electrons_photonIso_);
      tree_->SetBranchStatus("electrons.puIso", 1); tree_->SetBranchAddress("electrons.puIso", &electrons_puIso_);
      tree_->SetBranchStatus("electrons.dB", 1); tree_->SetBranchAddress("electrons.dB", &electrons_dB_);
      tree_->SetBranchStatus("electrons.ipDXY", 1); tree_->SetBranchAddress("electrons.ipDXY", &electrons_ipDXY_);
      tree_->SetBranchStatus("electrons.dz", 1); tree_->SetBranchAddress("electrons.dz", &electrons_dz_);
      tree_->SetBranchStatus("electrons.nMissingInnerHits", 1); tree_->SetBranchAddress("electrons.nMissingInnerHits", &electrons_nMissingInnerHits_);
      tree_->SetBranchStatus("electrons.r9", 1); tree_->SetBranchAddress("electrons.r9", &electrons_r9_);
      tree_->SetBranchStatus("electrons.ESCOverETrack", 1); tree_->SetBranchAddress("electrons.ESCOverETrack", &electrons_ESCOverETrack_);
      tree_->SetBranchStatus("electrons.DEtaSCTrk", 1); tree_->SetBranchAddress("electrons.DEtaSCTrk", &electrons_DEtaSCTrk_);
      tree_->SetBranchStatus("electrons.DPhiSCTrk", 1); tree_->SetBranchAddress("electrons.DPhiSCTrk", &electrons_DPhiSCTrk_);
      tree_->SetBranchStatus("electrons.ecalEnergy", 1); tree_->SetBranchAddress("electrons.ecalEnergy", &electrons_ecalEnergy_);
      tree_->SetBranchStatus("electrons.passConversionVeto", 1); tree_->SetBranchAddress("electrons.passConversionVeto", &electrons_passConversionVeto_);
      tree_->SetBranchStatus("electrons.isEB", 1); tree_->SetBranchAddress("electrons.isEB", &electrons_isEB_);
      tree_->SetBranchStatus("electrons.isEE", 1); tree_->SetBranchAddress("electrons.isEE", &electrons_isEE_);
      tree_->SetBranchStatus("electrons.isEBGap", 1); tree_->SetBranchAddress("electrons.isEBGap", &electrons_isEBGap_);
      tree_->SetBranchStatus("electrons.isEBEtaGap", 1); tree_->SetBranchAddress("electrons.isEBEtaGap", &electrons_isEBEtaGap_);
      tree_->SetBranchStatus("electrons.isEBPhiGap", 1); tree_->SetBranchAddress("electrons.isEBPhiGap", &electrons_isEBPhiGap_);
      tree_->SetBranchStatus("electrons.isEEGap", 1); tree_->SetBranchAddress("electrons.isEEGap", &electrons_isEEGap_);
      tree_->SetBranchStatus("electrons.isEERingGap", 1); tree_->SetBranchAddress("electrons.isEERingGap", &electrons_isEERingGap_);
      tree_->SetBranchStatus("electrons.isEEDeeGap", 1); tree_->SetBranchAddress("electrons.isEEDeeGap", &electrons_isEEDeeGap_);
      tree_->SetBranchStatus("electrons.isEBEEGap", 1); tree_->SetBranchAddress("electrons.isEBEEGap", &electrons_isEBEEGap_);
      tree_->SetBranchStatus("electrons.isElectron", 1); tree_->SetBranchAddress("electrons.isElectron", &electrons_isElectron_);
      tree_->SetBranchStatus("electrons.ecalSeed", 1); tree_->SetBranchAddress("electrons.ecalSeed", &electrons_ecalSeed_);
      tree_->SetBranchStatus("electrons.trackSeed", 1); tree_->SetBranchAddress("electrons.trackSeed", &electrons_trackSeed_);
      tree_->SetBranchStatus("electrons.eidCutLoose", 1); tree_->SetBranchAddress("electrons.eidCutLoose", &electrons_eidCutLoose_);
      tree_->SetBranchStatus("electrons.eidCutMedium", 1); tree_->SetBranchAddress("electrons.eidCutMedium", &electrons_eidCutMedium_);
      tree_->SetBranchStatus("electrons.eidCutTight", 1); tree_->SetBranchAddress("electrons.eidCutTight", &electrons_eidCutTight_);
      tree_->SetBranchStatus("electrons.eidCutVeto", 1); tree_->SetBranchAddress("electrons.eidCutVeto", &electrons_eidCutVeto_);
      tree_->SetBranchStatus("electrons.eidMVAWP80", 1); tree_->SetBranchAddress("electrons.eidMVAWP80", &electrons_eidMVAWP80_);
      tree_->SetBranchStatus("electrons.eidMVAWP90", 1); tree_->SetBranchAddress("electrons.eidMVAWP90", &electrons_eidMVAWP90_);
      tree_->SetBranchStatus("electrons.eidTrgMVAWP80", 1); tree_->SetBranchAddress("electrons.eidTrgMVAWP80", &electrons_eidTrgMVAWP80_);
      tree_->SetBranchStatus("electrons.eidTrgMVAWP90", 1); tree_->SetBranchAddress("electrons.eidTrgMVAWP90", &electrons_eidTrgMVAWP90_);
      tree_->SetBranchStatus("electrons.pfHadronIso", 1); tree_->SetBranchAddress("electrons.pfHadronIso", &electrons_pfHadronIso_);
      tree_->SetBranchStatus("electrons.pfNeutralIso", 1); tree_->SetBranchAddress("electrons.pfNeutralIso", &electrons_pfNeutralIso_);
      tree_->SetBranchStatus("electrons.pfPhotonIso", 1); tree_->SetBranchAddress("electrons.pfPhotonIso", &electrons_pfPhotonIso_);
      tree_->SetBranchStatus("electrons.HLT_Ele27_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("electrons.HLT_Ele27_eta2p1_WPLoose_Gsf", &electrons_HLT_Ele27_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("electrons.HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("electrons.HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &electrons_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("electrons.HLT_Ele22_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("electrons.HLT_Ele22_eta2p1_WPLoose_Gsf", &electrons_HLT_Ele22_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("electrons.HLT_Ele23_WPLoose_Gsf", 1); tree_->SetBranchAddress("electrons.HLT_Ele23_WPLoose_Gsf", &electrons_HLT_Ele23_WPLoose_Gsf_);
      tree_->SetBranchStatus("electrons.HLT_Ele27_WPLoose_Gsf", 1); tree_->SetBranchAddress("electrons.HLT_Ele27_WPLoose_Gsf", &electrons_HLT_Ele27_WPLoose_Gsf_);
      tree_->SetBranchStatus("electrons.e1x5", 1); tree_->SetBranchAddress("electrons.e1x5", &electrons_e1x5_);
      tree_->SetBranchStatus("electrons.e5x5", 1); tree_->SetBranchAddress("electrons.e5x5", &electrons_e5x5_);
      tree_->SetBranchStatus("electrons.sigmaIEtaIEta", 1); tree_->SetBranchAddress("electrons.sigmaIEtaIEta", &electrons_sigmaIEtaIEta_);
      tree_->SetBranchStatus("electrons.full5x5_sigmaIEtaIEta", 1); tree_->SetBranchAddress("electrons.full5x5_sigmaIEtaIEta", &electrons_full5x5_sigmaIEtaIEta_);
      tree_->SetBranchStatus("electrons.sigmaIPhiIPhi", 1); tree_->SetBranchAddress("electrons.sigmaIPhiIPhi", &electrons_sigmaIPhiIPhi_);
      tree_->SetBranchStatus("electrons.hadronicOverEM", 1); tree_->SetBranchAddress("electrons.hadronicOverEM", &electrons_hadronicOverEM_);
      tree_->SetBranchStatus("electrons.x", 1); tree_->SetBranchAddress("electrons.x", &electrons_x_);
      tree_->SetBranchStatus("electrons.y", 1); tree_->SetBranchAddress("electrons.y", &electrons_y_);
      tree_->SetBranchStatus("electrons.z", 1); tree_->SetBranchAddress("electrons.z", &electrons_z_);
      tree_->SetBranchStatus("electrons.energy", 1); tree_->SetBranchAddress("electrons.energy", &electrons_energy_);
      tree_->SetBranchStatus("electrons.rawEnergy", 1); tree_->SetBranchAddress("electrons.rawEnergy", &electrons_rawEnergy_);
      tree_->SetBranchStatus("electrons.phiWidth", 1); tree_->SetBranchAddress("electrons.phiWidth", &electrons_phiWidth_);
      tree_->SetBranchStatus("electrons.etaWidth", 1); tree_->SetBranchAddress("electrons.etaWidth", &electrons_etaWidth_);
      are_electrons_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenjets(){
    if(!are_genjets_loaded_){
      tree_->SetBranchStatus("genjets.pt", 1); tree_->SetBranchAddress("genjets.pt", &genjets_pt_);
      tree_->SetBranchStatus("genjets.eta", 1); tree_->SetBranchAddress("genjets.eta", &genjets_eta_);
      tree_->SetBranchStatus("genjets.phi", 1); tree_->SetBranchAddress("genjets.phi", &genjets_phi_);
      tree_->SetBranchStatus("genjets.charge", 1); tree_->SetBranchAddress("genjets.charge", &genjets_charge_);
      tree_->SetBranchStatus("genjets.e", 1); tree_->SetBranchAddress("genjets.e", &genjets_e_);
      tree_->SetBranchStatus("genjets.invisibleEnergy", 1); tree_->SetBranchAddress("genjets.invisibleEnergy", &genjets_invisibleEnergy_);
      tree_->SetBranchStatus("genjets.pdgId", 1); tree_->SetBranchAddress("genjets.pdgId", &genjets_pdgId_);
      are_genjets_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadRho(){
    if(!are_rho_loaded_){
      tree_->SetBranchStatus("rho.value", 1); tree_->SetBranchAddress("rho.value", &rho_value_);
      are_rho_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadJets(){
    if(!are_jets_loaded_){
      tree_->SetBranchStatus("jets.pt", 1); tree_->SetBranchAddress("jets.pt", &jets_pt_);
      tree_->SetBranchStatus("jets.eta", 1); tree_->SetBranchAddress("jets.eta", &jets_eta_);
      tree_->SetBranchStatus("jets.phi", 1); tree_->SetBranchAddress("jets.phi", &jets_phi_);
      tree_->SetBranchStatus("jets.charge", 1); tree_->SetBranchAddress("jets.charge", &jets_charge_);
      tree_->SetBranchStatus("jets.e", 1); tree_->SetBranchAddress("jets.e", &jets_e_);
      tree_->SetBranchStatus("jets.area", 1); tree_->SetBranchAddress("jets.area", &jets_area_);
      tree_->SetBranchStatus("jets.mass", 1); tree_->SetBranchAddress("jets.mass", &jets_mass_);
      tree_->SetBranchStatus("jets.JESUnc", 1); tree_->SetBranchAddress("jets.JESUnc", &jets_JESUnc_);
      tree_->SetBranchStatus("jets.JER", 1); tree_->SetBranchAddress("jets.JER", &jets_JER_);
      tree_->SetBranchStatus("jets.JERUp", 1); tree_->SetBranchAddress("jets.JERUp", &jets_JERUp_);
      tree_->SetBranchStatus("jets.JERDown", 1); tree_->SetBranchAddress("jets.JERDown", &jets_JERDown_);
      tree_->SetBranchStatus("jets.uncorrPt", 1); tree_->SetBranchAddress("jets.uncorrPt", &jets_uncorrPt_);
      tree_->SetBranchStatus("jets.uncorrEta", 1); tree_->SetBranchAddress("jets.uncorrEta", &jets_uncorrEta_);
      tree_->SetBranchStatus("jets.uncorrPhi", 1); tree_->SetBranchAddress("jets.uncorrPhi", &jets_uncorrPhi_);
      tree_->SetBranchStatus("jets.uncorrM", 1); tree_->SetBranchAddress("jets.uncorrM", &jets_uncorrM_);
      tree_->SetBranchStatus("jets.uncorrEnergy", 1); tree_->SetBranchAddress("jets.uncorrEnergy", &jets_uncorrEnergy_);
      tree_->SetBranchStatus("jets.chargedHadronEnergyFraction", 1); tree_->SetBranchAddress("jets.chargedHadronEnergyFraction", &jets_chargedHadronEnergyFraction_);
      tree_->SetBranchStatus("jets.neutralHadronEnergyFraction", 1); tree_->SetBranchAddress("jets.neutralHadronEnergyFraction", &jets_neutralHadronEnergyFraction_);
      tree_->SetBranchStatus("jets.chargedEmEnergyFraction", 1); tree_->SetBranchAddress("jets.chargedEmEnergyFraction", &jets_chargedEmEnergyFraction_);
      tree_->SetBranchStatus("jets.neutralEmEnergyFraction", 1); tree_->SetBranchAddress("jets.neutralEmEnergyFraction", &jets_neutralEmEnergyFraction_);
      tree_->SetBranchStatus("jets.HFHadronEnergyFraction", 1); tree_->SetBranchAddress("jets.HFHadronEnergyFraction", &jets_HFHadronEnergyFraction_);
      tree_->SetBranchStatus("jets.HFEMEnergyFraction", 1); tree_->SetBranchAddress("jets.HFEMEnergyFraction", &jets_HFEMEnergyFraction_);
      tree_->SetBranchStatus("jets.muonEnergyFraction", 1); tree_->SetBranchAddress("jets.muonEnergyFraction", &jets_muonEnergyFraction_);
      tree_->SetBranchStatus("jets.chargedMultiplicity", 1); tree_->SetBranchAddress("jets.chargedMultiplicity", &jets_chargedMultiplicity_);
      tree_->SetBranchStatus("jets.neutralMultiplicity", 1); tree_->SetBranchAddress("jets.neutralMultiplicity", &jets_neutralMultiplicity_);
      tree_->SetBranchStatus("jets.numChargedHadrons", 1); tree_->SetBranchAddress("jets.numChargedHadrons", &jets_numChargedHadrons_);
      tree_->SetBranchStatus("jets.numNeutralHadrons", 1); tree_->SetBranchAddress("jets.numNeutralHadrons", &jets_numNeutralHadrons_);
      tree_->SetBranchStatus("jets.numPhotons", 1); tree_->SetBranchAddress("jets.numPhotons", &jets_numPhotons_);
      tree_->SetBranchStatus("jets.numElectrons", 1); tree_->SetBranchAddress("jets.numElectrons", &jets_numElectrons_);
      tree_->SetBranchStatus("jets.numMuons", 1); tree_->SetBranchAddress("jets.numMuons", &jets_numMuons_);
      tree_->SetBranchStatus("jets.numForwardEMs", 1); tree_->SetBranchAddress("jets.numForwardEMs", &jets_numForwardEMs_);
      tree_->SetBranchStatus("jets.numForwardHads", 1); tree_->SetBranchAddress("jets.numForwardHads", &jets_numForwardHads_);
      tree_->SetBranchStatus("jets.numberOfDaughters", 1); tree_->SetBranchAddress("jets.numberOfDaughters", &jets_numberOfDaughters_);
      tree_->SetBranchStatus("jets.puId", 1); tree_->SetBranchAddress("jets.puId", &jets_puId_);
      tree_->SetBranchStatus("jets.jetBProb", 1); tree_->SetBranchAddress("jets.jetBProb", &jets_jetBProb_);
      tree_->SetBranchStatus("jets.jetProb", 1); tree_->SetBranchAddress("jets.jetProb", &jets_jetProb_);
      tree_->SetBranchStatus("jets.trkHiPur", 1); tree_->SetBranchAddress("jets.trkHiPur", &jets_trkHiPur_);
      tree_->SetBranchStatus("jets.trkHiEff", 1); tree_->SetBranchAddress("jets.trkHiEff", &jets_trkHiEff_);
      tree_->SetBranchStatus("jets.ssvHiEff", 1); tree_->SetBranchAddress("jets.ssvHiEff", &jets_ssvHiEff_);
      tree_->SetBranchStatus("jets.ssvHiPur", 1); tree_->SetBranchAddress("jets.ssvHiPur", &jets_ssvHiPur_);
      tree_->SetBranchStatus("jets.csv", 1); tree_->SetBranchAddress("jets.csv", &jets_csv_);
      tree_->SetBranchStatus("jets.csvIncl", 1); tree_->SetBranchAddress("jets.csvIncl", &jets_csvIncl_);
      tree_->SetBranchStatus("jets.CvsLtag", 1); tree_->SetBranchAddress("jets.CvsLtag", &jets_CvsLtag_);
      tree_->SetBranchStatus("jets.CombinedMVA", 1); tree_->SetBranchAddress("jets.CombinedMVA", &jets_CombinedMVA_);
      tree_->SetBranchStatus("jets.CvsBtag", 1); tree_->SetBranchAddress("jets.CvsBtag", &jets_CvsBtag_);
      tree_->SetBranchStatus("jets.vtxMass", 1); tree_->SetBranchAddress("jets.vtxMass", &jets_vtxMass_);
      tree_->SetBranchStatus("jets.vtxNtracks", 1); tree_->SetBranchAddress("jets.vtxNtracks", &jets_vtxNtracks_);
      tree_->SetBranchStatus("jets.vtx3DVal", 1); tree_->SetBranchAddress("jets.vtx3DVal", &jets_vtx3DVal_);
      tree_->SetBranchStatus("jets.vtx3DSig", 1); tree_->SetBranchAddress("jets.vtx3DSig", &jets_vtx3DSig_);
      tree_->SetBranchStatus("jets.partonFlavour", 1); tree_->SetBranchAddress("jets.partonFlavour", &jets_partonFlavour_);
      tree_->SetBranchStatus("jets.hadronFlavour", 1); tree_->SetBranchAddress("jets.hadronFlavour", &jets_hadronFlavour_);
      are_jets_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadMuons(){
    if(!are_muons_loaded_){
      tree_->SetBranchStatus("muons.pt", 1); tree_->SetBranchAddress("muons.pt", &muons_pt_);
      tree_->SetBranchStatus("muons.eta", 1); tree_->SetBranchAddress("muons.eta", &muons_eta_);
      tree_->SetBranchStatus("muons.phi", 1); tree_->SetBranchAddress("muons.phi", &muons_phi_);
      tree_->SetBranchStatus("muons.charge", 1); tree_->SetBranchAddress("muons.charge", &muons_charge_);
      tree_->SetBranchStatus("muons.dB", 1); tree_->SetBranchAddress("muons.dB", &muons_dB_);
      tree_->SetBranchStatus("muons.ipDXY", 1); tree_->SetBranchAddress("muons.ipDXY", &muons_ipDXY_);
      tree_->SetBranchStatus("muons.dz", 1); tree_->SetBranchAddress("muons.dz", &muons_dz_);
      tree_->SetBranchStatus("muons.nMissingInnerHits", 1); tree_->SetBranchAddress("muons.nMissingInnerHits", &muons_nMissingInnerHits_);
      tree_->SetBranchStatus("muons.chargedIso", 1); tree_->SetBranchAddress("muons.chargedIso", &muons_chargedIso_);
      tree_->SetBranchStatus("muons.neutralIso", 1); tree_->SetBranchAddress("muons.neutralIso", &muons_neutralIso_);
      tree_->SetBranchStatus("muons.photonIso", 1); tree_->SetBranchAddress("muons.photonIso", &muons_photonIso_);
      tree_->SetBranchStatus("muons.puIso", 1); tree_->SetBranchAddress("muons.puIso", &muons_puIso_);
      tree_->SetBranchStatus("muons.ECalEnergy", 1); tree_->SetBranchAddress("muons.ECalEnergy", &muons_ECalEnergy_);
      tree_->SetBranchStatus("muons.HCalEnergy", 1); tree_->SetBranchAddress("muons.HCalEnergy", &muons_HCalEnergy_);
      tree_->SetBranchStatus("muons.numChambers", 1); tree_->SetBranchAddress("muons.numChambers", &muons_numChambers_);
      tree_->SetBranchStatus("muons.numMatchedStations", 1); tree_->SetBranchAddress("muons.numMatchedStations", &muons_numMatchedStations_);
      tree_->SetBranchStatus("muons.trackiso", 1); tree_->SetBranchAddress("muons.trackiso", &muons_trackiso_);
      tree_->SetBranchStatus("muons.ecaliso", 1); tree_->SetBranchAddress("muons.ecaliso", &muons_ecaliso_);
      tree_->SetBranchStatus("muons.hcaliso", 1); tree_->SetBranchAddress("muons.hcaliso", &muons_hcaliso_);
      tree_->SetBranchStatus("muons.pfChargedIso04", 1); tree_->SetBranchAddress("muons.pfChargedIso04", &muons_pfChargedIso04_);
      tree_->SetBranchStatus("muons.pfNeutralIso04", 1); tree_->SetBranchAddress("muons.pfNeutralIso04", &muons_pfNeutralIso04_);
      tree_->SetBranchStatus("muons.pfPhotonIso04", 1); tree_->SetBranchAddress("muons.pfPhotonIso04", &muons_pfPhotonIso04_);
      tree_->SetBranchStatus("muons.pfPUIso04", 1); tree_->SetBranchAddress("muons.pfPUIso04", &muons_pfPUIso04_);
      tree_->SetBranchStatus("muons.trkIso03", 1); tree_->SetBranchAddress("muons.trkIso03", &muons_trkIso03_);
      tree_->SetBranchStatus("muons.ptErr", 1); tree_->SetBranchAddress("muons.ptErr", &muons_ptErr_);
      tree_->SetBranchStatus("muons.chi2", 1); tree_->SetBranchAddress("muons.chi2", &muons_chi2_);
      tree_->SetBranchStatus("muons.ndof", 1); tree_->SetBranchAddress("muons.ndof", &muons_ndof_);
      tree_->SetBranchStatus("muons.validHits", 1); tree_->SetBranchAddress("muons.validHits", &muons_validHits_);
      tree_->SetBranchStatus("muons.pixelHits", 1); tree_->SetBranchAddress("muons.pixelHits", &muons_pixelHits_);
      tree_->SetBranchStatus("muons.trackerLayers", 1); tree_->SetBranchAddress("muons.trackerLayers", &muons_trackerLayers_);
      tree_->SetBranchStatus("muons.isGlobal", 1); tree_->SetBranchAddress("muons.isGlobal", &muons_isGlobal_);
      tree_->SetBranchStatus("muons.isTracker", 1); tree_->SetBranchAddress("muons.isTracker", &muons_isTracker_);
      tree_->SetBranchStatus("muons.isCalo", 1); tree_->SetBranchAddress("muons.isCalo", &muons_isCalo_);
      tree_->SetBranchStatus("muons.isPF", 1); tree_->SetBranchAddress("muons.isPF", &muons_isPF_);
      tree_->SetBranchStatus("muons.isStandAlone", 1); tree_->SetBranchAddress("muons.isStandAlone", &muons_isStandAlone_);
      tree_->SetBranchStatus("muons.isLoose", 1); tree_->SetBranchAddress("muons.isLoose", &muons_isLoose_);
      tree_->SetBranchStatus("muons.HLT_DoubleIsoMu17_eta2p1", 1); tree_->SetBranchAddress("muons.HLT_DoubleIsoMu17_eta2p1", &muons_HLT_DoubleIsoMu17_eta2p1_);
      tree_->SetBranchStatus("muons.HLT_IsoMu24_eta2p1", 1); tree_->SetBranchAddress("muons.HLT_IsoMu24_eta2p1", &muons_HLT_IsoMu24_eta2p1_);
      tree_->SetBranchStatus("muons.HLT_IsoMu20_eta2p1", 1); tree_->SetBranchAddress("muons.HLT_IsoMu20_eta2p1", &muons_HLT_IsoMu20_eta2p1_);
      tree_->SetBranchStatus("muons.HLT_IsoMu18", 1); tree_->SetBranchAddress("muons.HLT_IsoMu18", &muons_HLT_IsoMu18_);
      tree_->SetBranchStatus("muons.HLT_IsoTkMu20", 1); tree_->SetBranchAddress("muons.HLT_IsoTkMu20", &muons_HLT_IsoTkMu20_);
      tree_->SetBranchStatus("muons.HLT_IsoMu20", 1); tree_->SetBranchAddress("muons.HLT_IsoMu20", &muons_HLT_IsoMu20_);
      tree_->SetBranchStatus("muons.HLT_IsoTkMu27", 1); tree_->SetBranchAddress("muons.HLT_IsoTkMu27", &muons_HLT_IsoTkMu27_);
      tree_->SetBranchStatus("muons.HLT_IsoMu27", 1); tree_->SetBranchAddress("muons.HLT_IsoMu27", &muons_HLT_IsoMu27_);
      tree_->SetBranchStatus("muons.HLT_Mu50", 1); tree_->SetBranchAddress("muons.HLT_Mu50", &muons_HLT_Mu50_);
      tree_->SetBranchStatus("muons.HLT_Mu45_eta2p1", 1); tree_->SetBranchAddress("muons.HLT_Mu45_eta2p1", &muons_HLT_Mu45_eta2p1_);
      tree_->SetBranchStatus("muons.HLT_IsoMu22", 1); tree_->SetBranchAddress("muons.HLT_IsoMu22", &muons_HLT_IsoMu22_);
      tree_->SetBranchStatus("muons.HLT_IsoTkMu22", 1); tree_->SetBranchAddress("muons.HLT_IsoTkMu22", &muons_HLT_IsoTkMu22_);
      are_muons_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenparticles(){
    if(!are_genParticles_loaded_){
      tree_->SetBranchStatus("genParticles.pt", 1); tree_->SetBranchAddress("genParticles.pt", &genParticles_pt_);
      tree_->SetBranchStatus("genParticles.eta", 1); tree_->SetBranchAddress("genParticles.eta", &genParticles_eta_);
      tree_->SetBranchStatus("genParticles.phi", 1); tree_->SetBranchAddress("genParticles.phi", &genParticles_phi_);
      tree_->SetBranchStatus("genParticles.charge", 1); tree_->SetBranchAddress("genParticles.charge", &genParticles_charge_);
      tree_->SetBranchStatus("genParticles.e", 1); tree_->SetBranchAddress("genParticles.e", &genParticles_e_);
      tree_->SetBranchStatus("genParticles.vx", 1); tree_->SetBranchAddress("genParticles.vx", &genParticles_vx_);
      tree_->SetBranchStatus("genParticles.vy", 1); tree_->SetBranchAddress("genParticles.vy", &genParticles_vy_);
      tree_->SetBranchStatus("genParticles.vz", 1); tree_->SetBranchAddress("genParticles.vz", &genParticles_vz_);
      tree_->SetBranchStatus("genParticles.pdgId", 1); tree_->SetBranchAddress("genParticles.pdgId", &genParticles_pdgId_);
      tree_->SetBranchStatus("genParticles.status", 1); tree_->SetBranchAddress("genParticles.status", &genParticles_status_);
      tree_->SetBranchStatus("genParticles.idx", 1); tree_->SetBranchAddress("genParticles.idx", &genParticles_idx_);
      tree_->SetBranchStatus("genParticles.momIdx", 1); tree_->SetBranchAddress("genParticles.momIdx", &genParticles_momIdx_);
      tree_->SetBranchStatus("genParticles.nDaught", 1); tree_->SetBranchAddress("genParticles.nDaught", &genParticles_nDaught_);
      tree_->SetBranchStatus("genParticles.firstDaughtIdx", 1); tree_->SetBranchAddress("genParticles.firstDaughtIdx", &genParticles_firstDaughtIdx_);
      are_genParticles_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadVertexs(){
    if(!are_vertexs_loaded_){
      tree_->SetBranchStatus("vertexs.x", 1); tree_->SetBranchAddress("vertexs.x", &vertexs_x_);
      tree_->SetBranchStatus("vertexs.y", 1); tree_->SetBranchAddress("vertexs.y", &vertexs_y_);
      tree_->SetBranchStatus("vertexs.z", 1); tree_->SetBranchAddress("vertexs.z", &vertexs_z_);
      tree_->SetBranchStatus("vertexs.chi2", 1); tree_->SetBranchAddress("vertexs.chi2", &vertexs_chi2_);
      tree_->SetBranchStatus("vertexs.ndof", 1); tree_->SetBranchAddress("vertexs.ndof", &vertexs_ndof_);
      tree_->SetBranchStatus("vertexs.nTracks", 1); tree_->SetBranchAddress("vertexs.nTracks", &vertexs_nTracks_);
      are_vertexs_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadPuinfos(){
    if(!are_PUInfos_loaded_){
      tree_->SetBranchStatus("PUInfos.bx", 1); tree_->SetBranchAddress("PUInfos.bx", &PUInfos_bx_);
      tree_->SetBranchStatus("PUInfos.nPU", 1); tree_->SetBranchAddress("PUInfos.nPU", &PUInfos_nPU_);
      tree_->SetBranchStatus("PUInfos.nInteractions", 1); tree_->SetBranchAddress("PUInfos.nInteractions", &PUInfos_nInteractions_);
      are_PUInfos_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  

  const Npnlolhe NPNLOLHE(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadNpnlolhe();
  
    Npnlolhe obj;
    obj.setnpnlo(NPNLOLHE_npnlo_);
  
    return obj;
  }
  
  const vector<Met>& METs(){
    if(METs_.size() > 0) return METs_;
    loadMets();
  	METs_.reserve(METs_px_->size());
    auto it_METs_px_ = METs_px_->cbegin();
    auto it_METs_py_ = METs_py_->cbegin();
    auto it_METs_pxsmear_ = METs_pxsmear_->cbegin();
    auto it_METs_pysmear_ = METs_pysmear_->cbegin();
    auto it_METs_pxunc_ = METs_pxunc_->cbegin();
    auto it_METs_pyunc_ = METs_pyunc_->cbegin();
    auto it_METs_pxuncJES_ = METs_pxuncJES_->cbegin();
    auto it_METs_pyuncJES_ = METs_pyuncJES_->cbegin();
    auto it_METs_pxuncJER_ = METs_pxuncJER_->cbegin();
    auto it_METs_pyuncJER_ = METs_pyuncJER_->cbegin();
    for(; it_METs_px_ != METs_px_->cend(); ){
      Met obj;
      obj.setpx(*it_METs_px_);
      obj.setpy(*it_METs_py_);
      obj.setpxsmear(*it_METs_pxsmear_);
      obj.setpysmear(*it_METs_pysmear_);
      obj.setpxunc(*it_METs_pxunc_);
      obj.setpyunc(*it_METs_pyunc_);
      obj.setpxuncJES(*it_METs_pxuncJES_);
      obj.setpyuncJES(*it_METs_pyuncJES_);
      obj.setpxuncJER(*it_METs_pxuncJER_);
      obj.setpyuncJER(*it_METs_pyuncJER_);
      
      METs_.push_back( obj );
      ++it_METs_px_;
      ++it_METs_py_;
      ++it_METs_pxsmear_;
      ++it_METs_pysmear_;
      ++it_METs_pxunc_;
      ++it_METs_pyunc_;
      ++it_METs_pxuncJES_;
      ++it_METs_pyuncJES_;
      ++it_METs_pxuncJER_;
      ++it_METs_pyuncJER_;
    }
    return METs_;
  }
  
  const vector<Mcweight>& MCWeights(){
    if(MCWeights_.size() > 0) return MCWeights_;
    loadMcweights();
  	MCWeights_.reserve(MCWeights_weights_->size());
    auto it_MCWeights_weights_ = MCWeights_weights_->cbegin();
    for(; it_MCWeights_weights_ != MCWeights_weights_->cend(); ){
      Mcweight obj;
      obj.setweights(*it_MCWeights_weights_);
      
      MCWeights_.push_back( obj );
      ++it_MCWeights_weights_;
    }
    return MCWeights_;
  }
  
  const Geninfo genInfo(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGeninfo();
  
    Geninfo obj;
    obj.setweight(genInfo_weight_);
    obj.setpdfid1(genInfo_pdfid1_);
    obj.setpdfid2(genInfo_pdfid2_);
    obj.setx1(genInfo_x1_);
    obj.setx2(genInfo_x2_);
    obj.setrenScale(genInfo_renScale_);
  
    return obj;
  }
  
  const vector<Photon>& photons(){
    if(photons_.size() > 0) return photons_;
    loadPhotons();
  	photons_.reserve(photons_pt_->size());
    auto it_photons_pt_ = photons_pt_->cbegin();
    auto it_photons_eta_ = photons_eta_->cbegin();
    auto it_photons_phi_ = photons_phi_->cbegin();
    auto it_photons_charge_ = photons_charge_->cbegin();
    auto it_photons_x_ = photons_x_->cbegin();
    auto it_photons_y_ = photons_y_->cbegin();
    auto it_photons_z_ = photons_z_->cbegin();
    auto it_photons_energy_ = photons_energy_->cbegin();
    auto it_photons_rawEnergy_ = photons_rawEnergy_->cbegin();
    auto it_photons_phiWidth_ = photons_phiWidth_->cbegin();
    auto it_photons_etaWidth_ = photons_etaWidth_->cbegin();
    auto it_photons_e3x3_ = photons_e3x3_->cbegin();
    auto it_photons_maxCrystalEnergy_ = photons_maxCrystalEnergy_->cbegin();
    auto it_photons_isEB_ = photons_isEB_->cbegin();
    auto it_photons_isEE_ = photons_isEE_->cbegin();
    auto it_photons_isPFlowPhoton_ = photons_isPFlowPhoton_->cbegin();
    auto it_photons_hasConversionTracks_ = photons_hasConversionTracks_->cbegin();
    auto it_photons_hasPixelSeed_ = photons_hasPixelSeed_->cbegin();
    for(; it_photons_pt_ != photons_pt_->cend(); ){
      Photon obj;
      obj.setcharge(*it_photons_charge_);
      obj.setx(*it_photons_x_);
      obj.sety(*it_photons_y_);
      obj.setz(*it_photons_z_);
      obj.setenergy(*it_photons_energy_);
      obj.setrawEnergy(*it_photons_rawEnergy_);
      obj.setphiWidth(*it_photons_phiWidth_);
      obj.setetaWidth(*it_photons_etaWidth_);
      obj.sete3x3(*it_photons_e3x3_);
      obj.setmaxCrystalEnergy(*it_photons_maxCrystalEnergy_);
      obj.setisEB(*it_photons_isEB_);
      obj.setisEE(*it_photons_isEE_);
      obj.setisPFlowPhoton(*it_photons_isPFlowPhoton_);
      obj.sethasConversionTracks(*it_photons_hasConversionTracks_);
      obj.sethasPixelSeed(*it_photons_hasPixelSeed_);
      obj.setLotentzVector(*it_photons_pt_, *it_photons_eta_, *it_photons_phi_);
      photons_.push_back( obj );
      ++it_photons_pt_;
      ++it_photons_eta_;
      ++it_photons_phi_;
      ++it_photons_charge_;
      ++it_photons_x_;
      ++it_photons_y_;
      ++it_photons_z_;
      ++it_photons_energy_;
      ++it_photons_rawEnergy_;
      ++it_photons_phiWidth_;
      ++it_photons_etaWidth_;
      ++it_photons_e3x3_;
      ++it_photons_maxCrystalEnergy_;
      ++it_photons_isEB_;
      ++it_photons_isEE_;
      ++it_photons_isPFlowPhoton_;
      ++it_photons_hasConversionTracks_;
      ++it_photons_hasPixelSeed_;
    }
    return photons_;
  }
  
  const vector<Lhepaticle>& LHEPaticles(){
    if(LHEPaticles_.size() > 0) return LHEPaticles_;
    loadLhepaticles();
  	LHEPaticles_.reserve(LHEPaticles_px_->size());
    auto it_LHEPaticles_px_ = LHEPaticles_px_->cbegin();
    auto it_LHEPaticles_py_ = LHEPaticles_py_->cbegin();
    auto it_LHEPaticles_pz_ = LHEPaticles_pz_->cbegin();
    auto it_LHEPaticles_e_ = LHEPaticles_e_->cbegin();
    auto it_LHEPaticles_pdgid_ = LHEPaticles_pdgid_->cbegin();
    auto it_LHEPaticles_status_ = LHEPaticles_status_->cbegin();
    auto it_LHEPaticles_fmother_ = LHEPaticles_fmother_->cbegin();
    auto it_LHEPaticles_lmother_ = LHEPaticles_lmother_->cbegin();
    for(; it_LHEPaticles_px_ != LHEPaticles_px_->cend(); ){
      Lhepaticle obj;
      obj.setpx(*it_LHEPaticles_px_);
      obj.setpy(*it_LHEPaticles_py_);
      obj.setpz(*it_LHEPaticles_pz_);
      obj.sete(*it_LHEPaticles_e_);
      obj.setpdgid(*it_LHEPaticles_pdgid_);
      obj.setstatus(*it_LHEPaticles_status_);
      obj.setfmother(*it_LHEPaticles_fmother_);
      obj.setlmother(*it_LHEPaticles_lmother_);
      
      LHEPaticles_.push_back( obj );
      ++it_LHEPaticles_px_;
      ++it_LHEPaticles_py_;
      ++it_LHEPaticles_pz_;
      ++it_LHEPaticles_e_;
      ++it_LHEPaticles_pdgid_;
      ++it_LHEPaticles_status_;
      ++it_LHEPaticles_fmother_;
      ++it_LHEPaticles_lmother_;
    }
    return LHEPaticles_;
  }
  
  const Filter filter(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadFilter();
  
    Filter obj;
    obj.setFlag_goodVertices(filter_Flag_goodVertices_);
    obj.setFlag_CSCTightHaloFilter(filter_Flag_CSCTightHaloFilter_);
    obj.setFlag_trkPOGFilters(filter_Flag_trkPOGFilters_);
    obj.setFlag_trkPOG_logErrorTooManyClusters(filter_Flag_trkPOG_logErrorTooManyClusters_);
    obj.setFlag_EcalDeadCellTriggerPrimitiveFilter(filter_Flag_EcalDeadCellTriggerPrimitiveFilter_);
    obj.setFlag_ecalLaserCorrFilter(filter_Flag_ecalLaserCorrFilter_);
    obj.setFlag_trkPOG_manystripclus53X(filter_Flag_trkPOG_manystripclus53X_);
    obj.setFlag_eeBadScFilter(filter_Flag_eeBadScFilter_);
    obj.setFlag_METFilters(filter_Flag_METFilters_);
    obj.setFlag_HBHENoiseFilter(filter_Flag_HBHENoiseFilter_);
    obj.setFlag_trkPOG_toomanystripclus53X(filter_Flag_trkPOG_toomanystripclus53X_);
    obj.setFlag_hcalLaserEventFilter(filter_Flag_hcalLaserEventFilter_);
    obj.setFlag_HBHENoiseIsoFilter(filter_Flag_HBHENoiseIsoFilter_);
    obj.setFlag_CSCTightHalo2015Filter(filter_Flag_CSCTightHalo2015Filter_);
  
    return obj;
  }
  
  const Trigger trigger(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadTrigger();
  
    Trigger obj;
    obj.setHLT_Ele27_eta2p1_WPLoose_Gsf(trigger_HLT_Ele27_eta2p1_WPLoose_Gsf_);
    obj.setHLT_DoubleIsoMu17_eta2p1(trigger_HLT_DoubleIsoMu17_eta2p1_);
    obj.setHLT_notexists(trigger_HLT_notexists_);
    obj.setHLT_IsoMu24_eta2p1(trigger_HLT_IsoMu24_eta2p1_);
    obj.setHLT_IsoMu20_eta2p1(trigger_HLT_IsoMu20_eta2p1_);
    obj.setHLT_DoubleEle33_CaloIdL_GsfTrkIdVL(trigger_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_);
    obj.setHLT_IsoMu18(trigger_HLT_IsoMu18_);
    obj.setHLT_IsoTkMu20(trigger_HLT_IsoTkMu20_);
    obj.setHLT_IsoMu20(trigger_HLT_IsoMu20_);
    obj.setHLT_IsoTkMu27(trigger_HLT_IsoTkMu27_);
    obj.setHLT_IsoMu27(trigger_HLT_IsoMu27_);
    obj.setHLT_Mu50(trigger_HLT_Mu50_);
    obj.setHLT_Mu45_eta2p1(trigger_HLT_Mu45_eta2p1_);
    obj.setHLT_Ele22_eta2p1_WPLoose_Gsf(trigger_HLT_Ele22_eta2p1_WPLoose_Gsf_);
    obj.setHLT_Ele23_WPLoose_Gsf(trigger_HLT_Ele23_WPLoose_Gsf_);
    obj.setHLT_Ele27_WPLoose_Gsf(trigger_HLT_Ele27_WPLoose_Gsf_);
    obj.setHLT_IsoMu22(trigger_HLT_IsoMu22_);
    obj.setHLT_IsoTkMu22(trigger_HLT_IsoTkMu22_);
  
    return obj;
  }
  
  const vector<Electron>& electrons(){
    if(electrons_.size() > 0) return electrons_;
    loadElectrons();
  	electrons_.reserve(electrons_pt_->size());
    auto it_electrons_pt_ = electrons_pt_->cbegin();
    auto it_electrons_eta_ = electrons_eta_->cbegin();
    auto it_electrons_phi_ = electrons_phi_->cbegin();
    auto it_electrons_charge_ = electrons_charge_->cbegin();
    auto it_electrons_chargedIso_ = electrons_chargedIso_->cbegin();
    auto it_electrons_neutralIso_ = electrons_neutralIso_->cbegin();
    auto it_electrons_photonIso_ = electrons_photonIso_->cbegin();
    auto it_electrons_puIso_ = electrons_puIso_->cbegin();
    auto it_electrons_dB_ = electrons_dB_->cbegin();
    auto it_electrons_ipDXY_ = electrons_ipDXY_->cbegin();
    auto it_electrons_dz_ = electrons_dz_->cbegin();
    auto it_electrons_nMissingInnerHits_ = electrons_nMissingInnerHits_->cbegin();
    auto it_electrons_r9_ = electrons_r9_->cbegin();
    auto it_electrons_ESCOverETrack_ = electrons_ESCOverETrack_->cbegin();
    auto it_electrons_DEtaSCTrk_ = electrons_DEtaSCTrk_->cbegin();
    auto it_electrons_DPhiSCTrk_ = electrons_DPhiSCTrk_->cbegin();
    auto it_electrons_ecalEnergy_ = electrons_ecalEnergy_->cbegin();
    auto it_electrons_passConversionVeto_ = electrons_passConversionVeto_->cbegin();
    auto it_electrons_isEB_ = electrons_isEB_->cbegin();
    auto it_electrons_isEE_ = electrons_isEE_->cbegin();
    auto it_electrons_isEBGap_ = electrons_isEBGap_->cbegin();
    auto it_electrons_isEBEtaGap_ = electrons_isEBEtaGap_->cbegin();
    auto it_electrons_isEBPhiGap_ = electrons_isEBPhiGap_->cbegin();
    auto it_electrons_isEEGap_ = electrons_isEEGap_->cbegin();
    auto it_electrons_isEERingGap_ = electrons_isEERingGap_->cbegin();
    auto it_electrons_isEEDeeGap_ = electrons_isEEDeeGap_->cbegin();
    auto it_electrons_isEBEEGap_ = electrons_isEBEEGap_->cbegin();
    auto it_electrons_isElectron_ = electrons_isElectron_->cbegin();
    auto it_electrons_ecalSeed_ = electrons_ecalSeed_->cbegin();
    auto it_electrons_trackSeed_ = electrons_trackSeed_->cbegin();
    auto it_electrons_eidCutLoose_ = electrons_eidCutLoose_->cbegin();
    auto it_electrons_eidCutMedium_ = electrons_eidCutMedium_->cbegin();
    auto it_electrons_eidCutTight_ = electrons_eidCutTight_->cbegin();
    auto it_electrons_eidCutVeto_ = electrons_eidCutVeto_->cbegin();
    auto it_electrons_eidMVAWP80_ = electrons_eidMVAWP80_->cbegin();
    auto it_electrons_eidMVAWP90_ = electrons_eidMVAWP90_->cbegin();
    auto it_electrons_eidTrgMVAWP80_ = electrons_eidTrgMVAWP80_->cbegin();
    auto it_electrons_eidTrgMVAWP90_ = electrons_eidTrgMVAWP90_->cbegin();
    auto it_electrons_pfHadronIso_ = electrons_pfHadronIso_->cbegin();
    auto it_electrons_pfNeutralIso_ = electrons_pfNeutralIso_->cbegin();
    auto it_electrons_pfPhotonIso_ = electrons_pfPhotonIso_->cbegin();
    auto it_electrons_HLT_Ele27_eta2p1_WPLoose_Gsf_ = electrons_HLT_Ele27_eta2p1_WPLoose_Gsf_->cbegin();
    auto it_electrons_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_ = electrons_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_->cbegin();
    auto it_electrons_HLT_Ele22_eta2p1_WPLoose_Gsf_ = electrons_HLT_Ele22_eta2p1_WPLoose_Gsf_->cbegin();
    auto it_electrons_HLT_Ele23_WPLoose_Gsf_ = electrons_HLT_Ele23_WPLoose_Gsf_->cbegin();
    auto it_electrons_HLT_Ele27_WPLoose_Gsf_ = electrons_HLT_Ele27_WPLoose_Gsf_->cbegin();
    auto it_electrons_e1x5_ = electrons_e1x5_->cbegin();
    auto it_electrons_e5x5_ = electrons_e5x5_->cbegin();
    auto it_electrons_sigmaIEtaIEta_ = electrons_sigmaIEtaIEta_->cbegin();
    auto it_electrons_full5x5_sigmaIEtaIEta_ = electrons_full5x5_sigmaIEtaIEta_->cbegin();
    auto it_electrons_sigmaIPhiIPhi_ = electrons_sigmaIPhiIPhi_->cbegin();
    auto it_electrons_hadronicOverEM_ = electrons_hadronicOverEM_->cbegin();
    auto it_electrons_x_ = electrons_x_->cbegin();
    auto it_electrons_y_ = electrons_y_->cbegin();
    auto it_electrons_z_ = electrons_z_->cbegin();
    auto it_electrons_energy_ = electrons_energy_->cbegin();
    auto it_electrons_rawEnergy_ = electrons_rawEnergy_->cbegin();
    auto it_electrons_phiWidth_ = electrons_phiWidth_->cbegin();
    auto it_electrons_etaWidth_ = electrons_etaWidth_->cbegin();
    for(; it_electrons_pt_ != electrons_pt_->cend(); ){
      Electron obj;
      obj.setcharge(*it_electrons_charge_);
      obj.setchargedIso(*it_electrons_chargedIso_);
      obj.setneutralIso(*it_electrons_neutralIso_);
      obj.setphotonIso(*it_electrons_photonIso_);
      obj.setpuIso(*it_electrons_puIso_);
      obj.setdB(*it_electrons_dB_);
      obj.setipDXY(*it_electrons_ipDXY_);
      obj.setdz(*it_electrons_dz_);
      obj.setnMissingInnerHits(*it_electrons_nMissingInnerHits_);
      obj.setr9(*it_electrons_r9_);
      obj.setESCOverETrack(*it_electrons_ESCOverETrack_);
      obj.setDEtaSCTrk(*it_electrons_DEtaSCTrk_);
      obj.setDPhiSCTrk(*it_electrons_DPhiSCTrk_);
      obj.setecalEnergy(*it_electrons_ecalEnergy_);
      obj.setpassConversionVeto(*it_electrons_passConversionVeto_);
      obj.setisEB(*it_electrons_isEB_);
      obj.setisEE(*it_electrons_isEE_);
      obj.setisEBGap(*it_electrons_isEBGap_);
      obj.setisEBEtaGap(*it_electrons_isEBEtaGap_);
      obj.setisEBPhiGap(*it_electrons_isEBPhiGap_);
      obj.setisEEGap(*it_electrons_isEEGap_);
      obj.setisEERingGap(*it_electrons_isEERingGap_);
      obj.setisEEDeeGap(*it_electrons_isEEDeeGap_);
      obj.setisEBEEGap(*it_electrons_isEBEEGap_);
      obj.setisElectron(*it_electrons_isElectron_);
      obj.setecalSeed(*it_electrons_ecalSeed_);
      obj.settrackSeed(*it_electrons_trackSeed_);
      obj.seteidCutLoose(*it_electrons_eidCutLoose_);
      obj.seteidCutMedium(*it_electrons_eidCutMedium_);
      obj.seteidCutTight(*it_electrons_eidCutTight_);
      obj.seteidCutVeto(*it_electrons_eidCutVeto_);
      obj.seteidMVAWP80(*it_electrons_eidMVAWP80_);
      obj.seteidMVAWP90(*it_electrons_eidMVAWP90_);
      obj.seteidTrgMVAWP80(*it_electrons_eidTrgMVAWP80_);
      obj.seteidTrgMVAWP90(*it_electrons_eidTrgMVAWP90_);
      obj.setpfHadronIso(*it_electrons_pfHadronIso_);
      obj.setpfNeutralIso(*it_electrons_pfNeutralIso_);
      obj.setpfPhotonIso(*it_electrons_pfPhotonIso_);
      obj.setHLT_Ele27_eta2p1_WPLoose_Gsf(*it_electrons_HLT_Ele27_eta2p1_WPLoose_Gsf_);
      obj.setHLT_DoubleEle33_CaloIdL_GsfTrkIdVL(*it_electrons_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_);
      obj.setHLT_Ele22_eta2p1_WPLoose_Gsf(*it_electrons_HLT_Ele22_eta2p1_WPLoose_Gsf_);
      obj.setHLT_Ele23_WPLoose_Gsf(*it_electrons_HLT_Ele23_WPLoose_Gsf_);
      obj.setHLT_Ele27_WPLoose_Gsf(*it_electrons_HLT_Ele27_WPLoose_Gsf_);
      obj.sete1x5(*it_electrons_e1x5_);
      obj.sete5x5(*it_electrons_e5x5_);
      obj.setsigmaIEtaIEta(*it_electrons_sigmaIEtaIEta_);
      obj.setfull5x5_sigmaIEtaIEta(*it_electrons_full5x5_sigmaIEtaIEta_);
      obj.setsigmaIPhiIPhi(*it_electrons_sigmaIPhiIPhi_);
      obj.sethadronicOverEM(*it_electrons_hadronicOverEM_);
      obj.setx(*it_electrons_x_);
      obj.sety(*it_electrons_y_);
      obj.setz(*it_electrons_z_);
      obj.setenergy(*it_electrons_energy_);
      obj.setrawEnergy(*it_electrons_rawEnergy_);
      obj.setphiWidth(*it_electrons_phiWidth_);
      obj.setetaWidth(*it_electrons_etaWidth_);
      obj.setLotentzVector(*it_electrons_pt_, *it_electrons_eta_, *it_electrons_phi_);
      electrons_.push_back( obj );
      ++it_electrons_pt_;
      ++it_electrons_eta_;
      ++it_electrons_phi_;
      ++it_electrons_charge_;
      ++it_electrons_chargedIso_;
      ++it_electrons_neutralIso_;
      ++it_electrons_photonIso_;
      ++it_electrons_puIso_;
      ++it_electrons_dB_;
      ++it_electrons_ipDXY_;
      ++it_electrons_dz_;
      ++it_electrons_nMissingInnerHits_;
      ++it_electrons_r9_;
      ++it_electrons_ESCOverETrack_;
      ++it_electrons_DEtaSCTrk_;
      ++it_electrons_DPhiSCTrk_;
      ++it_electrons_ecalEnergy_;
      ++it_electrons_passConversionVeto_;
      ++it_electrons_isEB_;
      ++it_electrons_isEE_;
      ++it_electrons_isEBGap_;
      ++it_electrons_isEBEtaGap_;
      ++it_electrons_isEBPhiGap_;
      ++it_electrons_isEEGap_;
      ++it_electrons_isEERingGap_;
      ++it_electrons_isEEDeeGap_;
      ++it_electrons_isEBEEGap_;
      ++it_electrons_isElectron_;
      ++it_electrons_ecalSeed_;
      ++it_electrons_trackSeed_;
      ++it_electrons_eidCutLoose_;
      ++it_electrons_eidCutMedium_;
      ++it_electrons_eidCutTight_;
      ++it_electrons_eidCutVeto_;
      ++it_electrons_eidMVAWP80_;
      ++it_electrons_eidMVAWP90_;
      ++it_electrons_eidTrgMVAWP80_;
      ++it_electrons_eidTrgMVAWP90_;
      ++it_electrons_pfHadronIso_;
      ++it_electrons_pfNeutralIso_;
      ++it_electrons_pfPhotonIso_;
      ++it_electrons_HLT_Ele27_eta2p1_WPLoose_Gsf_;
      ++it_electrons_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_;
      ++it_electrons_HLT_Ele22_eta2p1_WPLoose_Gsf_;
      ++it_electrons_HLT_Ele23_WPLoose_Gsf_;
      ++it_electrons_HLT_Ele27_WPLoose_Gsf_;
      ++it_electrons_e1x5_;
      ++it_electrons_e5x5_;
      ++it_electrons_sigmaIEtaIEta_;
      ++it_electrons_full5x5_sigmaIEtaIEta_;
      ++it_electrons_sigmaIPhiIPhi_;
      ++it_electrons_hadronicOverEM_;
      ++it_electrons_x_;
      ++it_electrons_y_;
      ++it_electrons_z_;
      ++it_electrons_energy_;
      ++it_electrons_rawEnergy_;
      ++it_electrons_phiWidth_;
      ++it_electrons_etaWidth_;
    }
    return electrons_;
  }
  
  const vector<Genjet>& genjets(){
    if(genjets_.size() > 0) return genjets_;
    loadGenjets();
  	genjets_.reserve(genjets_pt_->size());
    auto it_genjets_pt_ = genjets_pt_->cbegin();
    auto it_genjets_eta_ = genjets_eta_->cbegin();
    auto it_genjets_phi_ = genjets_phi_->cbegin();
    auto it_genjets_charge_ = genjets_charge_->cbegin();
    auto it_genjets_e_ = genjets_e_->cbegin();
    auto it_genjets_invisibleEnergy_ = genjets_invisibleEnergy_->cbegin();
    auto it_genjets_pdgId_ = genjets_pdgId_->cbegin();
    for(; it_genjets_pt_ != genjets_pt_->cend(); ){
      Genjet obj;
      obj.setcharge(*it_genjets_charge_);
      obj.sete(*it_genjets_e_);
      obj.setinvisibleEnergy(*it_genjets_invisibleEnergy_);
      obj.setpdgId(*it_genjets_pdgId_);
      obj.setLotentzVector(*it_genjets_pt_, *it_genjets_eta_, *it_genjets_phi_);
      genjets_.push_back( obj );
      ++it_genjets_pt_;
      ++it_genjets_eta_;
      ++it_genjets_phi_;
      ++it_genjets_charge_;
      ++it_genjets_e_;
      ++it_genjets_invisibleEnergy_;
      ++it_genjets_pdgId_;
    }
    return genjets_;
  }
  
  const Rho rho(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadRho();
  
    Rho obj;
    obj.setvalue(rho_value_);
  
    return obj;
  }
  
  const vector<Jet>& jets(){
    if(jets_.size() > 0) return jets_;
    loadJets();
  	jets_.reserve(jets_pt_->size());
    auto it_jets_pt_ = jets_pt_->cbegin();
    auto it_jets_eta_ = jets_eta_->cbegin();
    auto it_jets_phi_ = jets_phi_->cbegin();
    auto it_jets_charge_ = jets_charge_->cbegin();
    auto it_jets_e_ = jets_e_->cbegin();
    auto it_jets_area_ = jets_area_->cbegin();
    auto it_jets_mass_ = jets_mass_->cbegin();
    auto it_jets_JESUnc_ = jets_JESUnc_->cbegin();
    auto it_jets_JER_ = jets_JER_->cbegin();
    auto it_jets_JERUp_ = jets_JERUp_->cbegin();
    auto it_jets_JERDown_ = jets_JERDown_->cbegin();
    auto it_jets_uncorrPt_ = jets_uncorrPt_->cbegin();
    auto it_jets_uncorrEta_ = jets_uncorrEta_->cbegin();
    auto it_jets_uncorrPhi_ = jets_uncorrPhi_->cbegin();
    auto it_jets_uncorrM_ = jets_uncorrM_->cbegin();
    auto it_jets_uncorrEnergy_ = jets_uncorrEnergy_->cbegin();
    auto it_jets_chargedHadronEnergyFraction_ = jets_chargedHadronEnergyFraction_->cbegin();
    auto it_jets_neutralHadronEnergyFraction_ = jets_neutralHadronEnergyFraction_->cbegin();
    auto it_jets_chargedEmEnergyFraction_ = jets_chargedEmEnergyFraction_->cbegin();
    auto it_jets_neutralEmEnergyFraction_ = jets_neutralEmEnergyFraction_->cbegin();
    auto it_jets_HFHadronEnergyFraction_ = jets_HFHadronEnergyFraction_->cbegin();
    auto it_jets_HFEMEnergyFraction_ = jets_HFEMEnergyFraction_->cbegin();
    auto it_jets_muonEnergyFraction_ = jets_muonEnergyFraction_->cbegin();
    auto it_jets_chargedMultiplicity_ = jets_chargedMultiplicity_->cbegin();
    auto it_jets_neutralMultiplicity_ = jets_neutralMultiplicity_->cbegin();
    auto it_jets_numChargedHadrons_ = jets_numChargedHadrons_->cbegin();
    auto it_jets_numNeutralHadrons_ = jets_numNeutralHadrons_->cbegin();
    auto it_jets_numPhotons_ = jets_numPhotons_->cbegin();
    auto it_jets_numElectrons_ = jets_numElectrons_->cbegin();
    auto it_jets_numMuons_ = jets_numMuons_->cbegin();
    auto it_jets_numForwardEMs_ = jets_numForwardEMs_->cbegin();
    auto it_jets_numForwardHads_ = jets_numForwardHads_->cbegin();
    auto it_jets_numberOfDaughters_ = jets_numberOfDaughters_->cbegin();
    auto it_jets_puId_ = jets_puId_->cbegin();
    auto it_jets_jetBProb_ = jets_jetBProb_->cbegin();
    auto it_jets_jetProb_ = jets_jetProb_->cbegin();
    auto it_jets_trkHiPur_ = jets_trkHiPur_->cbegin();
    auto it_jets_trkHiEff_ = jets_trkHiEff_->cbegin();
    auto it_jets_ssvHiEff_ = jets_ssvHiEff_->cbegin();
    auto it_jets_ssvHiPur_ = jets_ssvHiPur_->cbegin();
    auto it_jets_csv_ = jets_csv_->cbegin();
    auto it_jets_csvIncl_ = jets_csvIncl_->cbegin();
    auto it_jets_CvsLtag_ = jets_CvsLtag_->cbegin();
    auto it_jets_CombinedMVA_ = jets_CombinedMVA_->cbegin();
    auto it_jets_CvsBtag_ = jets_CvsBtag_->cbegin();
    auto it_jets_vtxMass_ = jets_vtxMass_->cbegin();
    auto it_jets_vtxNtracks_ = jets_vtxNtracks_->cbegin();
    auto it_jets_vtx3DVal_ = jets_vtx3DVal_->cbegin();
    auto it_jets_vtx3DSig_ = jets_vtx3DSig_->cbegin();
    auto it_jets_partonFlavour_ = jets_partonFlavour_->cbegin();
    auto it_jets_hadronFlavour_ = jets_hadronFlavour_->cbegin();
    for(; it_jets_pt_ != jets_pt_->cend(); ){
      Jet obj;
      obj.setcharge(*it_jets_charge_);
      obj.sete(*it_jets_e_);
      obj.setarea(*it_jets_area_);
      obj.setJESUnc(*it_jets_JESUnc_);
      obj.setJER(*it_jets_JER_);
      obj.setJERUp(*it_jets_JERUp_);
      obj.setJERDown(*it_jets_JERDown_);
      obj.setuncorrPt(*it_jets_uncorrPt_);
      obj.setuncorrEta(*it_jets_uncorrEta_);
      obj.setuncorrPhi(*it_jets_uncorrPhi_);
      obj.setuncorrM(*it_jets_uncorrM_);
      obj.setuncorrEnergy(*it_jets_uncorrEnergy_);
      obj.setchargedHadronEnergyFraction(*it_jets_chargedHadronEnergyFraction_);
      obj.setneutralHadronEnergyFraction(*it_jets_neutralHadronEnergyFraction_);
      obj.setchargedEmEnergyFraction(*it_jets_chargedEmEnergyFraction_);
      obj.setneutralEmEnergyFraction(*it_jets_neutralEmEnergyFraction_);
      obj.setHFHadronEnergyFraction(*it_jets_HFHadronEnergyFraction_);
      obj.setHFEMEnergyFraction(*it_jets_HFEMEnergyFraction_);
      obj.setmuonEnergyFraction(*it_jets_muonEnergyFraction_);
      obj.setchargedMultiplicity(*it_jets_chargedMultiplicity_);
      obj.setneutralMultiplicity(*it_jets_neutralMultiplicity_);
      obj.setnumChargedHadrons(*it_jets_numChargedHadrons_);
      obj.setnumNeutralHadrons(*it_jets_numNeutralHadrons_);
      obj.setnumPhotons(*it_jets_numPhotons_);
      obj.setnumElectrons(*it_jets_numElectrons_);
      obj.setnumMuons(*it_jets_numMuons_);
      obj.setnumForwardEMs(*it_jets_numForwardEMs_);
      obj.setnumForwardHads(*it_jets_numForwardHads_);
      obj.setnumberOfDaughters(*it_jets_numberOfDaughters_);
      obj.setpuId(*it_jets_puId_);
      obj.setjetBProb(*it_jets_jetBProb_);
      obj.setjetProb(*it_jets_jetProb_);
      obj.settrkHiPur(*it_jets_trkHiPur_);
      obj.settrkHiEff(*it_jets_trkHiEff_);
      obj.setssvHiEff(*it_jets_ssvHiEff_);
      obj.setssvHiPur(*it_jets_ssvHiPur_);
      obj.setcsv(*it_jets_csv_);
      obj.setcsvIncl(*it_jets_csvIncl_);
      obj.setCvsLtag(*it_jets_CvsLtag_);
      obj.setCombinedMVA(*it_jets_CombinedMVA_);
      obj.setCvsBtag(*it_jets_CvsBtag_);
      obj.setvtxMass(*it_jets_vtxMass_);
      obj.setvtxNtracks(*it_jets_vtxNtracks_);
      obj.setvtx3DVal(*it_jets_vtx3DVal_);
      obj.setvtx3DSig(*it_jets_vtx3DSig_);
      obj.setpartonFlavour(*it_jets_partonFlavour_);
      obj.sethadronFlavour(*it_jets_hadronFlavour_);
      obj.setLotentzVector(*it_jets_pt_, *it_jets_eta_, *it_jets_phi_, *it_jets_mass_);
      jets_.push_back( obj );
      ++it_jets_pt_;
      ++it_jets_eta_;
      ++it_jets_phi_;
      ++it_jets_charge_;
      ++it_jets_e_;
      ++it_jets_area_;
      ++it_jets_mass_;
      ++it_jets_JESUnc_;
      ++it_jets_JER_;
      ++it_jets_JERUp_;
      ++it_jets_JERDown_;
      ++it_jets_uncorrPt_;
      ++it_jets_uncorrEta_;
      ++it_jets_uncorrPhi_;
      ++it_jets_uncorrM_;
      ++it_jets_uncorrEnergy_;
      ++it_jets_chargedHadronEnergyFraction_;
      ++it_jets_neutralHadronEnergyFraction_;
      ++it_jets_chargedEmEnergyFraction_;
      ++it_jets_neutralEmEnergyFraction_;
      ++it_jets_HFHadronEnergyFraction_;
      ++it_jets_HFEMEnergyFraction_;
      ++it_jets_muonEnergyFraction_;
      ++it_jets_chargedMultiplicity_;
      ++it_jets_neutralMultiplicity_;
      ++it_jets_numChargedHadrons_;
      ++it_jets_numNeutralHadrons_;
      ++it_jets_numPhotons_;
      ++it_jets_numElectrons_;
      ++it_jets_numMuons_;
      ++it_jets_numForwardEMs_;
      ++it_jets_numForwardHads_;
      ++it_jets_numberOfDaughters_;
      ++it_jets_puId_;
      ++it_jets_jetBProb_;
      ++it_jets_jetProb_;
      ++it_jets_trkHiPur_;
      ++it_jets_trkHiEff_;
      ++it_jets_ssvHiEff_;
      ++it_jets_ssvHiPur_;
      ++it_jets_csv_;
      ++it_jets_csvIncl_;
      ++it_jets_CvsLtag_;
      ++it_jets_CombinedMVA_;
      ++it_jets_CvsBtag_;
      ++it_jets_vtxMass_;
      ++it_jets_vtxNtracks_;
      ++it_jets_vtx3DVal_;
      ++it_jets_vtx3DSig_;
      ++it_jets_partonFlavour_;
      ++it_jets_hadronFlavour_;
    }
    return jets_;
  }
  
  const vector<Muon>& muons(){
    if(muons_.size() > 0) return muons_;
    loadMuons();
  	muons_.reserve(muons_pt_->size());
    auto it_muons_pt_ = muons_pt_->cbegin();
    auto it_muons_eta_ = muons_eta_->cbegin();
    auto it_muons_phi_ = muons_phi_->cbegin();
    auto it_muons_charge_ = muons_charge_->cbegin();
    auto it_muons_dB_ = muons_dB_->cbegin();
    auto it_muons_ipDXY_ = muons_ipDXY_->cbegin();
    auto it_muons_dz_ = muons_dz_->cbegin();
    auto it_muons_nMissingInnerHits_ = muons_nMissingInnerHits_->cbegin();
    auto it_muons_chargedIso_ = muons_chargedIso_->cbegin();
    auto it_muons_neutralIso_ = muons_neutralIso_->cbegin();
    auto it_muons_photonIso_ = muons_photonIso_->cbegin();
    auto it_muons_puIso_ = muons_puIso_->cbegin();
    auto it_muons_ECalEnergy_ = muons_ECalEnergy_->cbegin();
    auto it_muons_HCalEnergy_ = muons_HCalEnergy_->cbegin();
    auto it_muons_numChambers_ = muons_numChambers_->cbegin();
    auto it_muons_numMatchedStations_ = muons_numMatchedStations_->cbegin();
    auto it_muons_trackiso_ = muons_trackiso_->cbegin();
    auto it_muons_ecaliso_ = muons_ecaliso_->cbegin();
    auto it_muons_hcaliso_ = muons_hcaliso_->cbegin();
    auto it_muons_pfChargedIso04_ = muons_pfChargedIso04_->cbegin();
    auto it_muons_pfNeutralIso04_ = muons_pfNeutralIso04_->cbegin();
    auto it_muons_pfPhotonIso04_ = muons_pfPhotonIso04_->cbegin();
    auto it_muons_pfPUIso04_ = muons_pfPUIso04_->cbegin();
    auto it_muons_trkIso03_ = muons_trkIso03_->cbegin();
    auto it_muons_ptErr_ = muons_ptErr_->cbegin();
    auto it_muons_chi2_ = muons_chi2_->cbegin();
    auto it_muons_ndof_ = muons_ndof_->cbegin();
    auto it_muons_validHits_ = muons_validHits_->cbegin();
    auto it_muons_pixelHits_ = muons_pixelHits_->cbegin();
    auto it_muons_trackerLayers_ = muons_trackerLayers_->cbegin();
    auto it_muons_isGlobal_ = muons_isGlobal_->cbegin();
    auto it_muons_isTracker_ = muons_isTracker_->cbegin();
    auto it_muons_isCalo_ = muons_isCalo_->cbegin();
    auto it_muons_isPF_ = muons_isPF_->cbegin();
    auto it_muons_isStandAlone_ = muons_isStandAlone_->cbegin();
    auto it_muons_isLoose_ = muons_isLoose_->cbegin();
    auto it_muons_HLT_DoubleIsoMu17_eta2p1_ = muons_HLT_DoubleIsoMu17_eta2p1_->cbegin();
    auto it_muons_HLT_IsoMu24_eta2p1_ = muons_HLT_IsoMu24_eta2p1_->cbegin();
    auto it_muons_HLT_IsoMu20_eta2p1_ = muons_HLT_IsoMu20_eta2p1_->cbegin();
    auto it_muons_HLT_IsoMu18_ = muons_HLT_IsoMu18_->cbegin();
    auto it_muons_HLT_IsoTkMu20_ = muons_HLT_IsoTkMu20_->cbegin();
    auto it_muons_HLT_IsoMu20_ = muons_HLT_IsoMu20_->cbegin();
    auto it_muons_HLT_IsoTkMu27_ = muons_HLT_IsoTkMu27_->cbegin();
    auto it_muons_HLT_IsoMu27_ = muons_HLT_IsoMu27_->cbegin();
    auto it_muons_HLT_Mu50_ = muons_HLT_Mu50_->cbegin();
    auto it_muons_HLT_Mu45_eta2p1_ = muons_HLT_Mu45_eta2p1_->cbegin();
    auto it_muons_HLT_IsoMu22_ = muons_HLT_IsoMu22_->cbegin();
    auto it_muons_HLT_IsoTkMu22_ = muons_HLT_IsoTkMu22_->cbegin();
    for(; it_muons_pt_ != muons_pt_->cend(); ){
      Muon obj;
      obj.setcharge(*it_muons_charge_);
      obj.setdB(*it_muons_dB_);
      obj.setipDXY(*it_muons_ipDXY_);
      obj.setdz(*it_muons_dz_);
      obj.setnMissingInnerHits(*it_muons_nMissingInnerHits_);
      obj.setchargedIso(*it_muons_chargedIso_);
      obj.setneutralIso(*it_muons_neutralIso_);
      obj.setphotonIso(*it_muons_photonIso_);
      obj.setpuIso(*it_muons_puIso_);
      obj.setECalEnergy(*it_muons_ECalEnergy_);
      obj.setHCalEnergy(*it_muons_HCalEnergy_);
      obj.setnumChambers(*it_muons_numChambers_);
      obj.setnumMatchedStations(*it_muons_numMatchedStations_);
      obj.settrackiso(*it_muons_trackiso_);
      obj.setecaliso(*it_muons_ecaliso_);
      obj.sethcaliso(*it_muons_hcaliso_);
      obj.setpfChargedIso04(*it_muons_pfChargedIso04_);
      obj.setpfNeutralIso04(*it_muons_pfNeutralIso04_);
      obj.setpfPhotonIso04(*it_muons_pfPhotonIso04_);
      obj.setpfPUIso04(*it_muons_pfPUIso04_);
      obj.settrkIso03(*it_muons_trkIso03_);
      obj.setptErr(*it_muons_ptErr_);
      obj.setchi2(*it_muons_chi2_);
      obj.setndof(*it_muons_ndof_);
      obj.setvalidHits(*it_muons_validHits_);
      obj.setpixelHits(*it_muons_pixelHits_);
      obj.settrackerLayers(*it_muons_trackerLayers_);
      obj.setisGlobal(*it_muons_isGlobal_);
      obj.setisTracker(*it_muons_isTracker_);
      obj.setisCalo(*it_muons_isCalo_);
      obj.setisPF(*it_muons_isPF_);
      obj.setisStandAlone(*it_muons_isStandAlone_);
      obj.setisLoose(*it_muons_isLoose_);
      obj.setHLT_DoubleIsoMu17_eta2p1(*it_muons_HLT_DoubleIsoMu17_eta2p1_);
      obj.setHLT_IsoMu24_eta2p1(*it_muons_HLT_IsoMu24_eta2p1_);
      obj.setHLT_IsoMu20_eta2p1(*it_muons_HLT_IsoMu20_eta2p1_);
      obj.setHLT_IsoMu18(*it_muons_HLT_IsoMu18_);
      obj.setHLT_IsoTkMu20(*it_muons_HLT_IsoTkMu20_);
      obj.setHLT_IsoMu20(*it_muons_HLT_IsoMu20_);
      obj.setHLT_IsoTkMu27(*it_muons_HLT_IsoTkMu27_);
      obj.setHLT_IsoMu27(*it_muons_HLT_IsoMu27_);
      obj.setHLT_Mu50(*it_muons_HLT_Mu50_);
      obj.setHLT_Mu45_eta2p1(*it_muons_HLT_Mu45_eta2p1_);
      obj.setHLT_IsoMu22(*it_muons_HLT_IsoMu22_);
      obj.setHLT_IsoTkMu22(*it_muons_HLT_IsoTkMu22_);
      obj.setLotentzVector(*it_muons_pt_, *it_muons_eta_, *it_muons_phi_);
      muons_.push_back( obj );
      ++it_muons_pt_;
      ++it_muons_eta_;
      ++it_muons_phi_;
      ++it_muons_charge_;
      ++it_muons_dB_;
      ++it_muons_ipDXY_;
      ++it_muons_dz_;
      ++it_muons_nMissingInnerHits_;
      ++it_muons_chargedIso_;
      ++it_muons_neutralIso_;
      ++it_muons_photonIso_;
      ++it_muons_puIso_;
      ++it_muons_ECalEnergy_;
      ++it_muons_HCalEnergy_;
      ++it_muons_numChambers_;
      ++it_muons_numMatchedStations_;
      ++it_muons_trackiso_;
      ++it_muons_ecaliso_;
      ++it_muons_hcaliso_;
      ++it_muons_pfChargedIso04_;
      ++it_muons_pfNeutralIso04_;
      ++it_muons_pfPhotonIso04_;
      ++it_muons_pfPUIso04_;
      ++it_muons_trkIso03_;
      ++it_muons_ptErr_;
      ++it_muons_chi2_;
      ++it_muons_ndof_;
      ++it_muons_validHits_;
      ++it_muons_pixelHits_;
      ++it_muons_trackerLayers_;
      ++it_muons_isGlobal_;
      ++it_muons_isTracker_;
      ++it_muons_isCalo_;
      ++it_muons_isPF_;
      ++it_muons_isStandAlone_;
      ++it_muons_isLoose_;
      ++it_muons_HLT_DoubleIsoMu17_eta2p1_;
      ++it_muons_HLT_IsoMu24_eta2p1_;
      ++it_muons_HLT_IsoMu20_eta2p1_;
      ++it_muons_HLT_IsoMu18_;
      ++it_muons_HLT_IsoTkMu20_;
      ++it_muons_HLT_IsoMu20_;
      ++it_muons_HLT_IsoTkMu27_;
      ++it_muons_HLT_IsoMu27_;
      ++it_muons_HLT_Mu50_;
      ++it_muons_HLT_Mu45_eta2p1_;
      ++it_muons_HLT_IsoMu22_;
      ++it_muons_HLT_IsoTkMu22_;
    }
    return muons_;
  }
  
  const vector<Genparticle>& genParticles(){
    if(genParticles_.size() > 0) return genParticles_;
    loadGenparticles();
  	genParticles_.reserve(genParticles_pt_->size());
    auto it_genParticles_pt_ = genParticles_pt_->cbegin();
    auto it_genParticles_eta_ = genParticles_eta_->cbegin();
    auto it_genParticles_phi_ = genParticles_phi_->cbegin();
    auto it_genParticles_charge_ = genParticles_charge_->cbegin();
    auto it_genParticles_e_ = genParticles_e_->cbegin();
    auto it_genParticles_vx_ = genParticles_vx_->cbegin();
    auto it_genParticles_vy_ = genParticles_vy_->cbegin();
    auto it_genParticles_vz_ = genParticles_vz_->cbegin();
    auto it_genParticles_pdgId_ = genParticles_pdgId_->cbegin();
    auto it_genParticles_status_ = genParticles_status_->cbegin();
    auto it_genParticles_idx_ = genParticles_idx_->cbegin();
    auto it_genParticles_momIdx_ = genParticles_momIdx_->cbegin();
    auto it_genParticles_nDaught_ = genParticles_nDaught_->cbegin();
    auto it_genParticles_firstDaughtIdx_ = genParticles_firstDaughtIdx_->cbegin();
    for(; it_genParticles_pt_ != genParticles_pt_->cend(); ){
      Genparticle obj;
      obj.setcharge(*it_genParticles_charge_);
      obj.sete(*it_genParticles_e_);
      obj.setvx(*it_genParticles_vx_);
      obj.setvy(*it_genParticles_vy_);
      obj.setvz(*it_genParticles_vz_);
      obj.setpdgId(*it_genParticles_pdgId_);
      obj.setstatus(*it_genParticles_status_);
      obj.setidx(*it_genParticles_idx_);
      obj.setmomIdx(*it_genParticles_momIdx_);
      obj.setnDaught(*it_genParticles_nDaught_);
      obj.setfirstDaughtIdx(*it_genParticles_firstDaughtIdx_);
      obj.setLotentzVector(*it_genParticles_pt_, *it_genParticles_eta_, *it_genParticles_phi_);
      genParticles_.push_back( obj );
      ++it_genParticles_pt_;
      ++it_genParticles_eta_;
      ++it_genParticles_phi_;
      ++it_genParticles_charge_;
      ++it_genParticles_e_;
      ++it_genParticles_vx_;
      ++it_genParticles_vy_;
      ++it_genParticles_vz_;
      ++it_genParticles_pdgId_;
      ++it_genParticles_status_;
      ++it_genParticles_idx_;
      ++it_genParticles_momIdx_;
      ++it_genParticles_nDaught_;
      ++it_genParticles_firstDaughtIdx_;
    }
    return genParticles_;
  }
  
  const vector<Vertex>& vertexs(){
    if(vertexs_.size() > 0) return vertexs_;
    loadVertexs();
  	vertexs_.reserve(vertexs_x_->size());
    auto it_vertexs_x_ = vertexs_x_->cbegin();
    auto it_vertexs_y_ = vertexs_y_->cbegin();
    auto it_vertexs_z_ = vertexs_z_->cbegin();
    auto it_vertexs_chi2_ = vertexs_chi2_->cbegin();
    auto it_vertexs_ndof_ = vertexs_ndof_->cbegin();
    auto it_vertexs_nTracks_ = vertexs_nTracks_->cbegin();
    for(; it_vertexs_x_ != vertexs_x_->cend(); ){
      Vertex obj;
      obj.setx(*it_vertexs_x_);
      obj.sety(*it_vertexs_y_);
      obj.setz(*it_vertexs_z_);
      obj.setchi2(*it_vertexs_chi2_);
      obj.setndof(*it_vertexs_ndof_);
      obj.setnTracks(*it_vertexs_nTracks_);
      
      vertexs_.push_back( obj );
      ++it_vertexs_x_;
      ++it_vertexs_y_;
      ++it_vertexs_z_;
      ++it_vertexs_chi2_;
      ++it_vertexs_ndof_;
      ++it_vertexs_nTracks_;
    }
    return vertexs_;
  }
  
  const vector<Puinfo>& PUInfos(){
    if(PUInfos_.size() > 0) return PUInfos_;
    loadPuinfos();
  	PUInfos_.reserve(PUInfos_bx_->size());
    auto it_PUInfos_bx_ = PUInfos_bx_->cbegin();
    auto it_PUInfos_nPU_ = PUInfos_nPU_->cbegin();
    auto it_PUInfos_nInteractions_ = PUInfos_nInteractions_->cbegin();
    for(; it_PUInfos_bx_ != PUInfos_bx_->cend(); ){
      Puinfo obj;
      obj.setbx(*it_PUInfos_bx_);
      obj.setnPU(*it_PUInfos_nPU_);
      obj.setnInteractions(*it_PUInfos_nInteractions_);
      
      PUInfos_.push_back( obj );
      ++it_PUInfos_bx_;
      ++it_PUInfos_nPU_;
      ++it_PUInfos_nInteractions_;
    }
    return PUInfos_;
  }
  

private:
  TTree *tree_;
  Long64_t entries_;
  Long64_t current_entry_;
  Int_t trigger_HLT_Ele27_eta2p1_WPLoose_Gsf_;
  Int_t trigger_HLT_DoubleIsoMu17_eta2p1_;
  Int_t trigger_HLT_notexists_;
  Int_t trigger_HLT_IsoMu24_eta2p1_;
  Int_t trigger_HLT_IsoMu20_eta2p1_;
  Int_t trigger_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_;
  Int_t trigger_HLT_IsoMu18_;
  Int_t trigger_HLT_IsoTkMu20_;
  Int_t trigger_HLT_IsoMu20_;
  Int_t trigger_HLT_IsoTkMu27_;
  Int_t trigger_HLT_IsoMu27_;
  Int_t trigger_HLT_Mu50_;
  Int_t trigger_HLT_Mu45_eta2p1_;
  Int_t trigger_HLT_Ele22_eta2p1_WPLoose_Gsf_;
  Int_t trigger_HLT_Ele23_WPLoose_Gsf_;
  Int_t trigger_HLT_Ele27_WPLoose_Gsf_;
  Int_t trigger_HLT_IsoMu22_;
  Int_t trigger_HLT_IsoTkMu22_;
  Int_t filter_Flag_goodVertices_;
  Int_t filter_Flag_CSCTightHaloFilter_;
  Int_t filter_Flag_trkPOGFilters_;
  Int_t filter_Flag_trkPOG_logErrorTooManyClusters_;
  Int_t filter_Flag_EcalDeadCellTriggerPrimitiveFilter_;
  Int_t filter_Flag_ecalLaserCorrFilter_;
  Int_t filter_Flag_trkPOG_manystripclus53X_;
  Int_t filter_Flag_eeBadScFilter_;
  Int_t filter_Flag_METFilters_;
  Int_t filter_Flag_HBHENoiseFilter_;
  Int_t filter_Flag_trkPOG_toomanystripclus53X_;
  Int_t filter_Flag_hcalLaserEventFilter_;
  Int_t filter_Flag_HBHENoiseIsoFilter_;
  Int_t filter_Flag_CSCTightHalo2015Filter_;
  Double_t rho_value_;
  vector<float> *muons_pt_;
  vector<float> *muons_eta_;
  vector<float> *muons_phi_;
  vector<int> *muons_charge_;
  vector<float> *muons_dB_;
  vector<float> *muons_ipDXY_;
  vector<float> *muons_dz_;
  vector<float> *muons_nMissingInnerHits_;
  vector<float> *muons_chargedIso_;
  vector<float> *muons_neutralIso_;
  vector<float> *muons_photonIso_;
  vector<float> *muons_puIso_;
  vector<float> *muons_ECalEnergy_;
  vector<float> *muons_HCalEnergy_;
  vector<int> *muons_numChambers_;
  vector<int> *muons_numMatchedStations_;
  vector<float> *muons_trackiso_;
  vector<float> *muons_ecaliso_;
  vector<float> *muons_hcaliso_;
  vector<float> *muons_pfChargedIso04_;
  vector<float> *muons_pfNeutralIso04_;
  vector<float> *muons_pfPhotonIso04_;
  vector<float> *muons_pfPUIso04_;
  vector<float> *muons_trkIso03_;
  vector<float> *muons_ptErr_;
  vector<float> *muons_chi2_;
  vector<int> *muons_ndof_;
  vector<float> *muons_validHits_;
  vector<float> *muons_pixelHits_;
  vector<float> *muons_trackerLayers_;
  vector<bool> *muons_isGlobal_;
  vector<bool> *muons_isTracker_;
  vector<bool> *muons_isCalo_;
  vector<bool> *muons_isPF_;
  vector<bool> *muons_isStandAlone_;
  vector<bool> *muons_isLoose_;
  vector<bool> *muons_HLT_DoubleIsoMu17_eta2p1_;
  vector<bool> *muons_HLT_IsoMu24_eta2p1_;
  vector<bool> *muons_HLT_IsoMu20_eta2p1_;
  vector<bool> *muons_HLT_IsoMu18_;
  vector<bool> *muons_HLT_IsoTkMu20_;
  vector<bool> *muons_HLT_IsoMu20_;
  vector<bool> *muons_HLT_IsoTkMu27_;
  vector<bool> *muons_HLT_IsoMu27_;
  vector<bool> *muons_HLT_Mu50_;
  vector<bool> *muons_HLT_Mu45_eta2p1_;
  vector<bool> *muons_HLT_IsoMu22_;
  vector<bool> *muons_HLT_IsoTkMu22_;
  vector<float> *jets_pt_;
  vector<float> *jets_eta_;
  vector<float> *jets_phi_;
  vector<int> *jets_charge_;
  vector<float> *jets_e_;
  vector<float> *jets_area_;
  vector<float> *jets_mass_;
  vector<float> *jets_JESUnc_;
  vector<float> *jets_JER_;
  vector<float> *jets_JERUp_;
  vector<float> *jets_JERDown_;
  vector<float> *jets_uncorrPt_;
  vector<float> *jets_uncorrEta_;
  vector<float> *jets_uncorrPhi_;
  vector<float> *jets_uncorrM_;
  vector<float> *jets_uncorrEnergy_;
  vector<float> *jets_chargedHadronEnergyFraction_;
  vector<float> *jets_neutralHadronEnergyFraction_;
  vector<float> *jets_chargedEmEnergyFraction_;
  vector<float> *jets_neutralEmEnergyFraction_;
  vector<float> *jets_HFHadronEnergyFraction_;
  vector<float> *jets_HFEMEnergyFraction_;
  vector<float> *jets_muonEnergyFraction_;
  vector<float> *jets_chargedMultiplicity_;
  vector<float> *jets_neutralMultiplicity_;
  vector<float> *jets_numChargedHadrons_;
  vector<float> *jets_numNeutralHadrons_;
  vector<float> *jets_numPhotons_;
  vector<float> *jets_numElectrons_;
  vector<float> *jets_numMuons_;
  vector<float> *jets_numForwardEMs_;
  vector<float> *jets_numForwardHads_;
  vector<float> *jets_numberOfDaughters_;
  vector<float> *jets_puId_;
  vector<float> *jets_jetBProb_;
  vector<float> *jets_jetProb_;
  vector<float> *jets_trkHiPur_;
  vector<float> *jets_trkHiEff_;
  vector<float> *jets_ssvHiEff_;
  vector<float> *jets_ssvHiPur_;
  vector<float> *jets_csv_;
  vector<float> *jets_csvIncl_;
  vector<float> *jets_CvsLtag_;
  vector<float> *jets_CombinedMVA_;
  vector<float> *jets_CvsBtag_;
  vector<float> *jets_vtxMass_;
  vector<float> *jets_vtxNtracks_;
  vector<float> *jets_vtx3DVal_;
  vector<float> *jets_vtx3DSig_;
  vector<int> *jets_partonFlavour_;
  vector<int> *jets_hadronFlavour_;
  vector<float> *electrons_pt_;
  vector<float> *electrons_eta_;
  vector<float> *electrons_phi_;
  vector<int> *electrons_charge_;
  vector<float> *electrons_chargedIso_;
  vector<float> *electrons_neutralIso_;
  vector<float> *electrons_photonIso_;
  vector<float> *electrons_puIso_;
  vector<float> *electrons_dB_;
  vector<float> *electrons_ipDXY_;
  vector<float> *electrons_dz_;
  vector<float> *electrons_nMissingInnerHits_;
  vector<float> *electrons_r9_;
  vector<float> *electrons_ESCOverETrack_;
  vector<float> *electrons_DEtaSCTrk_;
  vector<float> *electrons_DPhiSCTrk_;
  vector<float> *electrons_ecalEnergy_;
  vector<bool> *electrons_passConversionVeto_;
  vector<bool> *electrons_isEB_;
  vector<bool> *electrons_isEE_;
  vector<bool> *electrons_isEBGap_;
  vector<bool> *electrons_isEBEtaGap_;
  vector<bool> *electrons_isEBPhiGap_;
  vector<bool> *electrons_isEEGap_;
  vector<bool> *electrons_isEERingGap_;
  vector<bool> *electrons_isEEDeeGap_;
  vector<bool> *electrons_isEBEEGap_;
  vector<bool> *electrons_isElectron_;
  vector<bool> *electrons_ecalSeed_;
  vector<bool> *electrons_trackSeed_;
  vector<float> *electrons_eidCutLoose_;
  vector<float> *electrons_eidCutMedium_;
  vector<float> *electrons_eidCutTight_;
  vector<float> *electrons_eidCutVeto_;
  vector<float> *electrons_eidMVAWP80_;
  vector<float> *electrons_eidMVAWP90_;
  vector<float> *electrons_eidTrgMVAWP80_;
  vector<float> *electrons_eidTrgMVAWP90_;
  vector<float> *electrons_pfHadronIso_;
  vector<float> *electrons_pfNeutralIso_;
  vector<float> *electrons_pfPhotonIso_;
  vector<bool> *electrons_HLT_Ele27_eta2p1_WPLoose_Gsf_;
  vector<bool> *electrons_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_;
  vector<bool> *electrons_HLT_Ele22_eta2p1_WPLoose_Gsf_;
  vector<bool> *electrons_HLT_Ele23_WPLoose_Gsf_;
  vector<bool> *electrons_HLT_Ele27_WPLoose_Gsf_;
  vector<float> *electrons_e1x5_;
  vector<float> *electrons_e5x5_;
  vector<float> *electrons_sigmaIEtaIEta_;
  vector<float> *electrons_full5x5_sigmaIEtaIEta_;
  vector<float> *electrons_sigmaIPhiIPhi_;
  vector<float> *electrons_hadronicOverEM_;
  vector<float> *electrons_x_;
  vector<float> *electrons_y_;
  vector<float> *electrons_z_;
  vector<float> *electrons_energy_;
  vector<float> *electrons_rawEnergy_;
  vector<float> *electrons_phiWidth_;
  vector<float> *electrons_etaWidth_;
  vector<float> *photons_pt_;
  vector<float> *photons_eta_;
  vector<float> *photons_phi_;
  vector<int> *photons_charge_;
  vector<float> *photons_x_;
  vector<float> *photons_y_;
  vector<float> *photons_z_;
  vector<float> *photons_energy_;
  vector<float> *photons_rawEnergy_;
  vector<float> *photons_phiWidth_;
  vector<float> *photons_etaWidth_;
  vector<float> *photons_e3x3_;
  vector<float> *photons_maxCrystalEnergy_;
  vector<bool> *photons_isEB_;
  vector<bool> *photons_isEE_;
  vector<bool> *photons_isPFlowPhoton_;
  vector<bool> *photons_hasConversionTracks_;
  vector<bool> *photons_hasPixelSeed_;
  vector<float> *vertexs_x_;
  vector<float> *vertexs_y_;
  vector<float> *vertexs_z_;
  vector<float> *vertexs_chi2_;
  vector<float> *vertexs_ndof_;
  vector<float> *vertexs_nTracks_;
  vector<float> *METs_px_;
  vector<float> *METs_py_;
  vector<float> *METs_pxsmear_;
  vector<float> *METs_pysmear_;
  vector<float> *METs_pxunc_;
  vector<float> *METs_pyunc_;
  vector<float> *METs_pxuncJES_;
  vector<float> *METs_pyuncJES_;
  vector<float> *METs_pxuncJER_;
  vector<float> *METs_pyuncJER_;
  Float_t genInfo_weight_;
  Float_t genInfo_pdfid1_;
  Float_t genInfo_pdfid2_;
  Float_t genInfo_x1_;
  Float_t genInfo_x2_;
  Float_t genInfo_renScale_;
  vector<float> *MCWeights_weights_;
  Int_t NPNLOLHE_npnlo_;
  vector<float> *PUInfos_bx_;
  vector<float> *PUInfos_nPU_;
  vector<float> *PUInfos_nInteractions_;
  vector<float> *genParticles_pt_;
  vector<float> *genParticles_eta_;
  vector<float> *genParticles_phi_;
  vector<int> *genParticles_charge_;
  vector<float> *genParticles_e_;
  vector<float> *genParticles_vx_;
  vector<float> *genParticles_vy_;
  vector<float> *genParticles_vz_;
  vector<int> *genParticles_pdgId_;
  vector<int> *genParticles_status_;
  vector<float> *LHEPaticles_px_;
  vector<float> *LHEPaticles_py_;
  vector<float> *LHEPaticles_pz_;
  vector<float> *LHEPaticles_e_;
  vector<int> *LHEPaticles_pdgid_;
  vector<int> *LHEPaticles_status_;
  vector<int> *LHEPaticles_fmother_;
  vector<int> *LHEPaticles_lmother_;
  vector<int> *genParticles_idx_;
  vector<vector<int> > *genParticles_momIdx_;
  vector<int> *genParticles_nDaught_;
  vector<int> *genParticles_firstDaughtIdx_;
  vector<float> *genjets_pt_;
  vector<float> *genjets_eta_;
  vector<float> *genjets_phi_;
  vector<int> *genjets_charge_;
  vector<float> *genjets_e_;
  vector<float> *genjets_invisibleEnergy_;
  vector<float> *genjets_pdgId_;
  bool are_NPNLOLHE_loaded_;
  Npnlolhe NPNLOLHE_;
  bool are_METs_loaded_;
  vector<Met> METs_;
  bool are_MCWeights_loaded_;
  vector<Mcweight> MCWeights_;
  bool are_genInfo_loaded_;
  Geninfo genInfo_;
  bool are_photons_loaded_;
  vector<Photon> photons_;
  bool are_LHEPaticles_loaded_;
  vector<Lhepaticle> LHEPaticles_;
  bool are_filter_loaded_;
  Filter filter_;
  bool are_trigger_loaded_;
  Trigger trigger_;
  bool are_electrons_loaded_;
  vector<Electron> electrons_;
  bool are_genjets_loaded_;
  vector<Genjet> genjets_;
  bool are_rho_loaded_;
  Rho rho_;
  bool are_jets_loaded_;
  vector<Jet> jets_;
  bool are_muons_loaded_;
  vector<Muon> muons_;
  bool are_genParticles_loaded_;
  vector<Genparticle> genParticles_;
  bool are_vertexs_loaded_;
  vector<Vertex> vertexs_;
  bool are_PUInfos_loaded_;
  vector<Puinfo> PUInfos_;
};

/*#include <iostream>
int test()
{
  TFile *f = TFile::Open("test_ntuple.root");
  TTree *t = (TTree*) f->Get("ntuple/events");
  URStreamer s(t);
  for(int i =0; i < 30; i++){
    s.next();
    vector<Muon> muons = s.muons();
    std::cout<< muons.size() << std::endl;
    for(int j=0; j<muons.size(); ++j){
      std::cout<< *(muons[j].pt) << "  ";
    }
    std::cout<< std::endl;
  }
  return 0;
  }*/

/*#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class vector+;
#pragma link C++ class URStreamer+;
{ LINK_OBJECTS }

#endif*/
#endif

