#ifndef URStreamer_h
#define URStreamer_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include <vector>
using namespace std;

class Electron: public TLorentzVector{
friend class URStreamer;
public:
//  Electron(const Float_t &i_deltaEtaSC_,const Float_t &i_dr03EcalRecHitSumEt_,const Float_t &i_dr03HcalDepth1TowerSumEt_,const Float_t &i_dr03TkSumPt_,const Float_t &i_dxy_,const Float_t &i_dxyErr_,const Float_t &i_dz_,const Float_t &i_dzErr_,const Float_t &i_eCorr_,const Float_t &i_eInvMinusPInv_,const Float_t &i_energyErr_,const Float_t &i_hoe_,const Float_t &i_ip3d_,const Float_t &i_miniPFRelIso_,const Float_t &i_miniPFRelIso_,const Float_t &i_mvaSpring16GP_,const Float_t &i_mvaSpring16HZZ_,const Float_t &i_pfRelIso03_,const Float_t &i_pfRelIso03_,const Float_t &i_r9_,const Float_t &i_sieie_,const Float_t &i_sip3d_,const Float_t &i_mvaTTH_,const Int_t &i_charge_,const Int_t &i_cutBased_,const Int_t &i_cutBased_,const Int_t &i_jetIdx_,const Int_t &i_pdgId_,const Int_t &i_photonIdx_,const Int_t &i_tightCharge_,const Int_t &i_vidNestedWPBitmap_,const Bool_t &i_convVeto_,const Bool_t &i_cutBased_,const Bool_t &i_isPFcand_,const UChar_t &i_lostHits_,const Bool_t &i_mvaSpring16GP_,const Bool_t &i_mvaSpring16GP_,const Bool_t &i_mvaSpring16HZZ_,const Int_t &i_genPartIdx_,const UChar_t &i_genPartFlav_,const UChar_t &i_cleanmask_):
//    
//  {}
  Electron():
    TLorentzVector(),
    deltaEtaSC_(0),
    dr03EcalRecHitSumEt_(0),
    dr03HcalDepth1TowerSumEt_(0),
    dr03TkSumPt_(0),
    dxy_(0),
    dxyErr_(0),
    dz_(0),
    dzErr_(0),
    eCorr_(0),
    eInvMinusPInv_(0),
    energyErr_(0),
    hoe_(0),
    ip3d_(0),
    miniPFRelIso_(0),
    miniPFRelIso_(0),
    mvaSpring16GP_(0),
    mvaSpring16HZZ_(0),
    pfRelIso03_(0),
    pfRelIso03_(0),
    r9_(0),
    sieie_(0),
    sip3d_(0),
    mvaTTH_(0),
    charge_(0),
    cutBased_(0),
    cutBased_(0),
    jetIdx_(0),
    pdgId_(0),
    photonIdx_(0),
    tightCharge_(0),
    vidNestedWPBitmap_(0),
    convVeto_(0),
    cutBased_(0),
    isPFcand_(0),
    lostHits_(0),
    mvaSpring16GP_(0),
    mvaSpring16GP_(0),
    mvaSpring16HZZ_(0),
    genPartIdx_(0),
    genPartFlav_(0),
    cleanmask_(0)
  {}
  Float_t deltaEtaSC() const {return deltaEtaSC_;}
  Float_t dr03EcalRecHitSumEt() const {return dr03EcalRecHitSumEt_;}
  Float_t dr03HcalDepth1TowerSumEt() const {return dr03HcalDepth1TowerSumEt_;}
  Float_t dr03TkSumPt() const {return dr03TkSumPt_;}
  Float_t dxy() const {return dxy_;}
  Float_t dxyErr() const {return dxyErr_;}
  Float_t dz() const {return dz_;}
  Float_t dzErr() const {return dzErr_;}
  Float_t eCorr() const {return eCorr_;}
  Float_t eInvMinusPInv() const {return eInvMinusPInv_;}
  Float_t energyErr() const {return energyErr_;}
  Float_t hoe() const {return hoe_;}
  Float_t ip3d() const {return ip3d_;}
  Float_t miniPFRelIso() const {return miniPFRelIso_;}
  Float_t miniPFRelIso() const {return miniPFRelIso_;}
  Float_t mvaSpring16GP() const {return mvaSpring16GP_;}
  Float_t mvaSpring16HZZ() const {return mvaSpring16HZZ_;}
  Float_t pfRelIso03() const {return pfRelIso03_;}
  Float_t pfRelIso03() const {return pfRelIso03_;}
  Float_t r9() const {return r9_;}
  Float_t sieie() const {return sieie_;}
  Float_t sip3d() const {return sip3d_;}
  Float_t mvaTTH() const {return mvaTTH_;}
  Int_t charge() const {return charge_;}
  Int_t cutBased() const {return cutBased_;}
  Int_t cutBased() const {return cutBased_;}
  Int_t jetIdx() const {return jetIdx_;}
  Int_t pdgId() const {return pdgId_;}
  Int_t photonIdx() const {return photonIdx_;}
  Int_t tightCharge() const {return tightCharge_;}
  Int_t vidNestedWPBitmap() const {return vidNestedWPBitmap_;}
  Bool_t convVeto() const {return convVeto_;}
  Bool_t cutBased() const {return cutBased_;}
  Bool_t isPFcand() const {return isPFcand_;}
  UChar_t lostHits() const {return lostHits_;}
  Bool_t mvaSpring16GP() const {return mvaSpring16GP_;}
  Bool_t mvaSpring16GP() const {return mvaSpring16GP_;}
  Bool_t mvaSpring16HZZ() const {return mvaSpring16HZZ_;}
  Int_t genPartIdx() const {return genPartIdx_;}
  UChar_t genPartFlav() const {return genPartFlav_;}
  UChar_t cleanmask() const {return cleanmask_;}
private:
  Float_t deltaEtaSC_;
  Float_t dr03EcalRecHitSumEt_;
  Float_t dr03HcalDepth1TowerSumEt_;
  Float_t dr03TkSumPt_;
  Float_t dxy_;
  Float_t dxyErr_;
  Float_t dz_;
  Float_t dzErr_;
  Float_t eCorr_;
  Float_t eInvMinusPInv_;
  Float_t energyErr_;
  Float_t hoe_;
  Float_t ip3d_;
  Float_t miniPFRelIso_;
  Float_t miniPFRelIso_;
  Float_t mvaSpring16GP_;
  Float_t mvaSpring16HZZ_;
  Float_t pfRelIso03_;
  Float_t pfRelIso03_;
  Float_t r9_;
  Float_t sieie_;
  Float_t sip3d_;
  Float_t mvaTTH_;
  Int_t charge_;
  Int_t cutBased_;
  Int_t cutBased_;
  Int_t jetIdx_;
  Int_t pdgId_;
  Int_t photonIdx_;
  Int_t tightCharge_;
  Int_t vidNestedWPBitmap_;
  Bool_t convVeto_;
  Bool_t cutBased_;
  Bool_t isPFcand_;
  UChar_t lostHits_;
  Bool_t mvaSpring16GP_;
  Bool_t mvaSpring16GP_;
  Bool_t mvaSpring16HZZ_;
  Int_t genPartIdx_;
  UChar_t genPartFlav_;
  UChar_t cleanmask_;
  void setdeltaEtaSC(const Float_t value) {deltaEtaSC_ = value;}
  void setdr03EcalRecHitSumEt(const Float_t value) {dr03EcalRecHitSumEt_ = value;}
  void setdr03HcalDepth1TowerSumEt(const Float_t value) {dr03HcalDepth1TowerSumEt_ = value;}
  void setdr03TkSumPt(const Float_t value) {dr03TkSumPt_ = value;}
  void setdxy(const Float_t value) {dxy_ = value;}
  void setdxyErr(const Float_t value) {dxyErr_ = value;}
  void setdz(const Float_t value) {dz_ = value;}
  void setdzErr(const Float_t value) {dzErr_ = value;}
  void seteCorr(const Float_t value) {eCorr_ = value;}
  void seteInvMinusPInv(const Float_t value) {eInvMinusPInv_ = value;}
  void setenergyErr(const Float_t value) {energyErr_ = value;}
  void sethoe(const Float_t value) {hoe_ = value;}
  void setip3d(const Float_t value) {ip3d_ = value;}
  void setminiPFRelIso(const Float_t value) {miniPFRelIso_ = value;}
  void setminiPFRelIso(const Float_t value) {miniPFRelIso_ = value;}
  void setmvaSpring16GP(const Float_t value) {mvaSpring16GP_ = value;}
  void setmvaSpring16HZZ(const Float_t value) {mvaSpring16HZZ_ = value;}
  void setpfRelIso03(const Float_t value) {pfRelIso03_ = value;}
  void setpfRelIso03(const Float_t value) {pfRelIso03_ = value;}
  void setr9(const Float_t value) {r9_ = value;}
  void setsieie(const Float_t value) {sieie_ = value;}
  void setsip3d(const Float_t value) {sip3d_ = value;}
  void setmvaTTH(const Float_t value) {mvaTTH_ = value;}
  void setcharge(const Int_t value) {charge_ = value;}
  void setcutBased(const Int_t value) {cutBased_ = value;}
  void setcutBased(const Int_t value) {cutBased_ = value;}
  void setjetIdx(const Int_t value) {jetIdx_ = value;}
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setphotonIdx(const Int_t value) {photonIdx_ = value;}
  void settightCharge(const Int_t value) {tightCharge_ = value;}
  void setvidNestedWPBitmap(const Int_t value) {vidNestedWPBitmap_ = value;}
  void setconvVeto(const Bool_t value) {convVeto_ = value;}
  void setcutBased(const Bool_t value) {cutBased_ = value;}
  void setisPFcand(const Bool_t value) {isPFcand_ = value;}
  void setlostHits(const UChar_t value) {lostHits_ = value;}
  void setmvaSpring16GP(const Bool_t value) {mvaSpring16GP_ = value;}
  void setmvaSpring16GP(const Bool_t value) {mvaSpring16GP_ = value;}
  void setmvaSpring16HZZ(const Bool_t value) {mvaSpring16HZZ_ = value;}
  void setgenPartIdx(const Int_t value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_t value) {genPartFlav_ = value;}
  void setcleanmask(const UChar_t value) {cleanmask_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Fatjet: public TLorentzVector{
friend class URStreamer;
public:
//  Fatjet(const Float_t &i_area_,const Float_t &i_btagCMVA_,const Float_t &i_btagCSVV2_,const Float_t &i_btagDeepB_,const Float_t &i_btagHbb_,const Float_t &i_msoftdrop_,const Float_t &i_msoftdrop_,const Float_t &i_n2b1_,const Float_t &i_n3b1_,const Float_t &i_tau1_,const Float_t &i_tau2_,const Float_t &i_tau3_,const Float_t &i_tau4_,const Int_t &i_jetId_,const Int_t &i_subJetIdx1_,const Int_t &i_subJetIdx2_):
//    
//  {}
  Fatjet():
    TLorentzVector(),
    area_(0),
    btagCMVA_(0),
    btagCSVV2_(0),
    btagDeepB_(0),
    btagHbb_(0),
    msoftdrop_(0),
    msoftdrop_(0),
    n2b1_(0),
    n3b1_(0),
    tau1_(0),
    tau2_(0),
    tau3_(0),
    tau4_(0),
    jetId_(0),
    subJetIdx1_(0),
    subJetIdx2_(0)
  {}
  Float_t area() const {return area_;}
  Float_t btagCMVA() const {return btagCMVA_;}
  Float_t btagCSVV2() const {return btagCSVV2_;}
  Float_t btagDeepB() const {return btagDeepB_;}
  Float_t btagHbb() const {return btagHbb_;}
  Float_t msoftdrop() const {return msoftdrop_;}
  Float_t msoftdrop() const {return msoftdrop_;}
  Float_t n2b1() const {return n2b1_;}
  Float_t n3b1() const {return n3b1_;}
  Float_t tau1() const {return tau1_;}
  Float_t tau2() const {return tau2_;}
  Float_t tau3() const {return tau3_;}
  Float_t tau4() const {return tau4_;}
  Int_t jetId() const {return jetId_;}
  Int_t subJetIdx1() const {return subJetIdx1_;}
  Int_t subJetIdx2() const {return subJetIdx2_;}
private:
  Float_t area_;
  Float_t btagCMVA_;
  Float_t btagCSVV2_;
  Float_t btagDeepB_;
  Float_t btagHbb_;
  Float_t msoftdrop_;
  Float_t msoftdrop_;
  Float_t n2b1_;
  Float_t n3b1_;
  Float_t tau1_;
  Float_t tau2_;
  Float_t tau3_;
  Float_t tau4_;
  Int_t jetId_;
  Int_t subJetIdx1_;
  Int_t subJetIdx2_;
  void setarea(const Float_t value) {area_ = value;}
  void setbtagCMVA(const Float_t value) {btagCMVA_ = value;}
  void setbtagCSVV2(const Float_t value) {btagCSVV2_ = value;}
  void setbtagDeepB(const Float_t value) {btagDeepB_ = value;}
  void setbtagHbb(const Float_t value) {btagHbb_ = value;}
  void setmsoftdrop(const Float_t value) {msoftdrop_ = value;}
  void setmsoftdrop(const Float_t value) {msoftdrop_ = value;}
  void setn2b1(const Float_t value) {n2b1_ = value;}
  void setn3b1(const Float_t value) {n3b1_ = value;}
  void settau1(const Float_t value) {tau1_ = value;}
  void settau2(const Float_t value) {tau2_ = value;}
  void settau3(const Float_t value) {tau3_ = value;}
  void settau4(const Float_t value) {tau4_ = value;}
  void setjetId(const Int_t value) {jetId_ = value;}
  void setsubJetIdx1(const Int_t value) {subJetIdx1_ = value;}
  void setsubJetIdx2(const Int_t value) {subJetIdx2_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Genjetak8: public TLorentzVector{
friend class URStreamer;
public:
//  Genjetak8(const Int_t &i_partonFlavour_,const UChar_t &i_hadronFlavour_):
//    
//  {}
  Genjetak8():
    TLorentzVector(),
    partonFlavour_(0),
    hadronFlavour_(0)
  {}
  Int_t partonFlavour() const {return partonFlavour_;}
  UChar_t hadronFlavour() const {return hadronFlavour_;}
private:
  Int_t partonFlavour_;
  UChar_t hadronFlavour_;
  void setpartonFlavour(const Int_t value) {partonFlavour_ = value;}
  void sethadronFlavour(const UChar_t value) {hadronFlavour_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Genjet: public TLorentzVector{
friend class URStreamer;
public:
//  Genjet(const Int_t &i_partonFlavour_,const UChar_t &i_hadronFlavour_,const Int_t &i_partonFlavour_,const UChar_t &i_hadronFlavour_):
//    
//  {}
  Genjet():
    TLorentzVector(),
    partonFlavour_(0),
    hadronFlavour_(0),
    partonFlavour_(0),
    hadronFlavour_(0)
  {}
  Int_t partonFlavour() const {return partonFlavour_;}
  UChar_t hadronFlavour() const {return hadronFlavour_;}
  Int_t partonFlavour() const {return partonFlavour_;}
  UChar_t hadronFlavour() const {return hadronFlavour_;}
private:
  Int_t partonFlavour_;
  UChar_t hadronFlavour_;
  Int_t partonFlavour_;
  UChar_t hadronFlavour_;
  void setpartonFlavour(const Int_t value) {partonFlavour_ = value;}
  void sethadronFlavour(const UChar_t value) {hadronFlavour_ = value;}
  void setpartonFlavour(const Int_t value) {partonFlavour_ = value;}
  void sethadronFlavour(const UChar_t value) {hadronFlavour_ = value;}
};

class Genpart: public TLorentzVector{
friend class URStreamer;
public:
//  Genpart(const Int_t &i_genPartIdxMother_,const Int_t &i_pdgId_,const Int_t &i_status_,const Int_t &i_statusFlags_):
//    
//  {}
  Genpart():
    TLorentzVector(),
    genPartIdxMother_(0),
    pdgId_(0),
    status_(0),
    statusFlags_(0)
  {}
  Int_t genPartIdxMother() const {return genPartIdxMother_;}
  Int_t pdgId() const {return pdgId_;}
  Int_t status() const {return status_;}
  Int_t statusFlags() const {return statusFlags_;}
private:
  Int_t genPartIdxMother_;
  Int_t pdgId_;
  Int_t status_;
  Int_t statusFlags_;
  void setgenPartIdxMother(const Int_t value) {genPartIdxMother_ = value;}
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setstatus(const Int_t value) {status_ = value;}
  void setstatusFlags(const Int_t value) {statusFlags_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Genvistau: public TLorentzVector{
friend class URStreamer;
public:
//  Genvistau(const Int_t &i_charge_,const Int_t &i_genPartIdxMother_,const Int_t &i_status_):
//    
//  {}
  Genvistau():
    TLorentzVector(),
    charge_(0),
    genPartIdxMother_(0),
    status_(0)
  {}
  Int_t charge() const {return charge_;}
  Int_t genPartIdxMother() const {return genPartIdxMother_;}
  Int_t status() const {return status_;}
private:
  Int_t charge_;
  Int_t genPartIdxMother_;
  Int_t status_;
  void setcharge(const Int_t value) {charge_ = value;}
  void setgenPartIdxMother(const Int_t value) {genPartIdxMother_ = value;}
  void setstatus(const Int_t value) {status_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Lhepdfweight{
friend class URStreamer;
public:
//  Lhepdfweight(const Float_t &i_LHEPdfWeight_):
//    
//  {}
  Lhepdfweight():
    LHEPdfWeight_(0)
  {}
  Float_t LHEPdfWeight() const {return LHEPdfWeight_;}
private:
  Float_t LHEPdfWeight_;
  void setLHEPdfWeight(const Float_t value) {LHEPdfWeight_ = value;}
};

class Lhescaleweight{
friend class URStreamer;
public:
//  Lhescaleweight(const Float_t &i_LHEScaleWeight_):
//    
//  {}
  Lhescaleweight():
    LHEScaleWeight_(0)
  {}
  Float_t LHEScaleWeight() const {return LHEScaleWeight_;}
private:
  Float_t LHEScaleWeight_;
  void setLHEScaleWeight(const Float_t value) {LHEScaleWeight_ = value;}
};

class Jet: public TLorentzVector{
friend class URStreamer;
public:
//  Jet(const Float_t &i_area_,const Float_t &i_btagCMVA_,const Float_t &i_btagCSVV2_,const Float_t &i_btagDeepB_,const Float_t &i_btagDeepC_,const Float_t &i_chEmEF_,const Float_t &i_chHEF_,const Float_t &i_neEmEF_,const Float_t &i_neHEF_,const Float_t &i_qgl_,const Float_t &i_rawFactor_,const Float_t &i_bReg_,const Int_t &i_electronIdx1_,const Int_t &i_electronIdx2_,const Int_t &i_jetId_,const Int_t &i_muonIdx1_,const Int_t &i_muonIdx2_,const Int_t &i_nConstituents_,const Int_t &i_nElectrons_,const Int_t &i_nMuons_,const Int_t &i_puId_,const Int_t &i_genJetIdx_,const Int_t &i_hadronFlavour_,const Int_t &i_partonFlavour_,const UChar_t &i_cleanmask_):
//    
//  {}
  Jet():
    TLorentzVector(),
    area_(0),
    btagCMVA_(0),
    btagCSVV2_(0),
    btagDeepB_(0),
    btagDeepC_(0),
    chEmEF_(0),
    chHEF_(0),
    neEmEF_(0),
    neHEF_(0),
    qgl_(0),
    rawFactor_(0),
    bReg_(0),
    electronIdx1_(0),
    electronIdx2_(0),
    jetId_(0),
    muonIdx1_(0),
    muonIdx2_(0),
    nConstituents_(0),
    nElectrons_(0),
    nMuons_(0),
    puId_(0),
    genJetIdx_(0),
    hadronFlavour_(0),
    partonFlavour_(0),
    cleanmask_(0)
  {}
  Float_t area() const {return area_;}
  Float_t btagCMVA() const {return btagCMVA_;}
  Float_t btagCSVV2() const {return btagCSVV2_;}
  Float_t btagDeepB() const {return btagDeepB_;}
  Float_t btagDeepC() const {return btagDeepC_;}
  Float_t chEmEF() const {return chEmEF_;}
  Float_t chHEF() const {return chHEF_;}
  Float_t neEmEF() const {return neEmEF_;}
  Float_t neHEF() const {return neHEF_;}
  Float_t qgl() const {return qgl_;}
  Float_t rawFactor() const {return rawFactor_;}
  Float_t bReg() const {return bReg_;}
  Int_t electronIdx1() const {return electronIdx1_;}
  Int_t electronIdx2() const {return electronIdx2_;}
  Int_t jetId() const {return jetId_;}
  Int_t muonIdx1() const {return muonIdx1_;}
  Int_t muonIdx2() const {return muonIdx2_;}
  Int_t nConstituents() const {return nConstituents_;}
  Int_t nElectrons() const {return nElectrons_;}
  Int_t nMuons() const {return nMuons_;}
  Int_t puId() const {return puId_;}
  Int_t genJetIdx() const {return genJetIdx_;}
  Int_t hadronFlavour() const {return hadronFlavour_;}
  Int_t partonFlavour() const {return partonFlavour_;}
  UChar_t cleanmask() const {return cleanmask_;}
private:
  Float_t area_;
  Float_t btagCMVA_;
  Float_t btagCSVV2_;
  Float_t btagDeepB_;
  Float_t btagDeepC_;
  Float_t chEmEF_;
  Float_t chHEF_;
  Float_t neEmEF_;
  Float_t neHEF_;
  Float_t qgl_;
  Float_t rawFactor_;
  Float_t bReg_;
  Int_t electronIdx1_;
  Int_t electronIdx2_;
  Int_t jetId_;
  Int_t muonIdx1_;
  Int_t muonIdx2_;
  Int_t nConstituents_;
  Int_t nElectrons_;
  Int_t nMuons_;
  Int_t puId_;
  Int_t genJetIdx_;
  Int_t hadronFlavour_;
  Int_t partonFlavour_;
  UChar_t cleanmask_;
  void setarea(const Float_t value) {area_ = value;}
  void setbtagCMVA(const Float_t value) {btagCMVA_ = value;}
  void setbtagCSVV2(const Float_t value) {btagCSVV2_ = value;}
  void setbtagDeepB(const Float_t value) {btagDeepB_ = value;}
  void setbtagDeepC(const Float_t value) {btagDeepC_ = value;}
  void setchEmEF(const Float_t value) {chEmEF_ = value;}
  void setchHEF(const Float_t value) {chHEF_ = value;}
  void setneEmEF(const Float_t value) {neEmEF_ = value;}
  void setneHEF(const Float_t value) {neHEF_ = value;}
  void setqgl(const Float_t value) {qgl_ = value;}
  void setrawFactor(const Float_t value) {rawFactor_ = value;}
  void setbReg(const Float_t value) {bReg_ = value;}
  void setelectronIdx1(const Int_t value) {electronIdx1_ = value;}
  void setelectronIdx2(const Int_t value) {electronIdx2_ = value;}
  void setjetId(const Int_t value) {jetId_ = value;}
  void setmuonIdx1(const Int_t value) {muonIdx1_ = value;}
  void setmuonIdx2(const Int_t value) {muonIdx2_ = value;}
  void setnConstituents(const Int_t value) {nConstituents_ = value;}
  void setnElectrons(const Int_t value) {nElectrons_ = value;}
  void setnMuons(const Int_t value) {nMuons_ = value;}
  void setpuId(const Int_t value) {puId_ = value;}
  void setgenJetIdx(const Int_t value) {genJetIdx_ = value;}
  void sethadronFlavour(const Int_t value) {hadronFlavour_ = value;}
  void setpartonFlavour(const Int_t value) {partonFlavour_ = value;}
  void setcleanmask(const UChar_t value) {cleanmask_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Muon: public TLorentzVector{
friend class URStreamer;
public:
//  Muon(const Float_t &i_dxy_,const Float_t &i_dxyErr_,const Float_t &i_dz_,const Float_t &i_dzErr_,const Float_t &i_ip3d_,const Float_t &i_miniPFRelIso_,const Float_t &i_miniPFRelIso_,const Float_t &i_pfRelIso03_,const Float_t &i_pfRelIso03_,const Float_t &i_pfRelIso04_,const Float_t &i_ptErr_,const Float_t &i_segmentComp_,const Float_t &i_sip3d_,const Float_t &i_mvaTTH_,const Int_t &i_charge_,const Int_t &i_jetIdx_,const Int_t &i_nStations_,const Int_t &i_nTrackerLayers_,const Int_t &i_pdgId_,const Int_t &i_tightCharge_,const UChar_t &i_highPtId_,const Bool_t &i_isPFcand_,const Bool_t &i_mediumId_,const Bool_t &i_softId_,const Bool_t &i_tightId_,const Int_t &i_genPartIdx_,const UChar_t &i_genPartFlav_,const UChar_t &i_cleanmask_):
//    
//  {}
  Muon():
    TLorentzVector(),
    dxy_(0),
    dxyErr_(0),
    dz_(0),
    dzErr_(0),
    ip3d_(0),
    miniPFRelIso_(0),
    miniPFRelIso_(0),
    pfRelIso03_(0),
    pfRelIso03_(0),
    pfRelIso04_(0),
    ptErr_(0),
    segmentComp_(0),
    sip3d_(0),
    mvaTTH_(0),
    charge_(0),
    jetIdx_(0),
    nStations_(0),
    nTrackerLayers_(0),
    pdgId_(0),
    tightCharge_(0),
    highPtId_(0),
    isPFcand_(0),
    mediumId_(0),
    softId_(0),
    tightId_(0),
    genPartIdx_(0),
    genPartFlav_(0),
    cleanmask_(0)
  {}
  Float_t dxy() const {return dxy_;}
  Float_t dxyErr() const {return dxyErr_;}
  Float_t dz() const {return dz_;}
  Float_t dzErr() const {return dzErr_;}
  Float_t ip3d() const {return ip3d_;}
  Float_t miniPFRelIso() const {return miniPFRelIso_;}
  Float_t miniPFRelIso() const {return miniPFRelIso_;}
  Float_t pfRelIso03() const {return pfRelIso03_;}
  Float_t pfRelIso03() const {return pfRelIso03_;}
  Float_t pfRelIso04() const {return pfRelIso04_;}
  Float_t ptErr() const {return ptErr_;}
  Float_t segmentComp() const {return segmentComp_;}
  Float_t sip3d() const {return sip3d_;}
  Float_t mvaTTH() const {return mvaTTH_;}
  Int_t charge() const {return charge_;}
  Int_t jetIdx() const {return jetIdx_;}
  Int_t nStations() const {return nStations_;}
  Int_t nTrackerLayers() const {return nTrackerLayers_;}
  Int_t pdgId() const {return pdgId_;}
  Int_t tightCharge() const {return tightCharge_;}
  UChar_t highPtId() const {return highPtId_;}
  Bool_t isPFcand() const {return isPFcand_;}
  Bool_t mediumId() const {return mediumId_;}
  Bool_t softId() const {return softId_;}
  Bool_t tightId() const {return tightId_;}
  Int_t genPartIdx() const {return genPartIdx_;}
  UChar_t genPartFlav() const {return genPartFlav_;}
  UChar_t cleanmask() const {return cleanmask_;}
private:
  Float_t dxy_;
  Float_t dxyErr_;
  Float_t dz_;
  Float_t dzErr_;
  Float_t ip3d_;
  Float_t miniPFRelIso_;
  Float_t miniPFRelIso_;
  Float_t pfRelIso03_;
  Float_t pfRelIso03_;
  Float_t pfRelIso04_;
  Float_t ptErr_;
  Float_t segmentComp_;
  Float_t sip3d_;
  Float_t mvaTTH_;
  Int_t charge_;
  Int_t jetIdx_;
  Int_t nStations_;
  Int_t nTrackerLayers_;
  Int_t pdgId_;
  Int_t tightCharge_;
  UChar_t highPtId_;
  Bool_t isPFcand_;
  Bool_t mediumId_;
  Bool_t softId_;
  Bool_t tightId_;
  Int_t genPartIdx_;
  UChar_t genPartFlav_;
  UChar_t cleanmask_;
  void setdxy(const Float_t value) {dxy_ = value;}
  void setdxyErr(const Float_t value) {dxyErr_ = value;}
  void setdz(const Float_t value) {dz_ = value;}
  void setdzErr(const Float_t value) {dzErr_ = value;}
  void setip3d(const Float_t value) {ip3d_ = value;}
  void setminiPFRelIso(const Float_t value) {miniPFRelIso_ = value;}
  void setminiPFRelIso(const Float_t value) {miniPFRelIso_ = value;}
  void setpfRelIso03(const Float_t value) {pfRelIso03_ = value;}
  void setpfRelIso03(const Float_t value) {pfRelIso03_ = value;}
  void setpfRelIso04(const Float_t value) {pfRelIso04_ = value;}
  void setptErr(const Float_t value) {ptErr_ = value;}
  void setsegmentComp(const Float_t value) {segmentComp_ = value;}
  void setsip3d(const Float_t value) {sip3d_ = value;}
  void setmvaTTH(const Float_t value) {mvaTTH_ = value;}
  void setcharge(const Int_t value) {charge_ = value;}
  void setjetIdx(const Int_t value) {jetIdx_ = value;}
  void setnStations(const Int_t value) {nStations_ = value;}
  void setnTrackerLayers(const Int_t value) {nTrackerLayers_ = value;}
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void settightCharge(const Int_t value) {tightCharge_ = value;}
  void sethighPtId(const UChar_t value) {highPtId_ = value;}
  void setisPFcand(const Bool_t value) {isPFcand_ = value;}
  void setmediumId(const Bool_t value) {mediumId_ = value;}
  void setsoftId(const Bool_t value) {softId_ = value;}
  void settightId(const Bool_t value) {tightId_ = value;}
  void setgenPartIdx(const Int_t value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_t value) {genPartFlav_ = value;}
  void setcleanmask(const UChar_t value) {cleanmask_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Photon: public TLorentzVector{
friend class URStreamer;
public:
//  Photon(const Float_t &i_eCorr_,const Float_t &i_energyErr_,const Float_t &i_hoe_,const Float_t &i_mvaID_,const Float_t &i_pfRelIso03_,const Float_t &i_pfRelIso03_,const Float_t &i_r9_,const Float_t &i_sieie_,const Int_t &i_charge_,const Int_t &i_cutBased_,const Int_t &i_electronIdx_,const Int_t &i_jetIdx_,const Int_t &i_pdgId_,const Int_t &i_vidNestedWPBitmap_,const Bool_t &i_electronVeto_,const Bool_t &i_mvaID_,const Bool_t &i_mvaID_,const Bool_t &i_pixelSeed_,const Int_t &i_genPartIdx_,const UChar_t &i_genPartFlav_,const UChar_t &i_cleanmask_):
//    
//  {}
  Photon():
    TLorentzVector(),
    eCorr_(0),
    energyErr_(0),
    hoe_(0),
    mvaID_(0),
    pfRelIso03_(0),
    pfRelIso03_(0),
    r9_(0),
    sieie_(0),
    charge_(0),
    cutBased_(0),
    electronIdx_(0),
    jetIdx_(0),
    pdgId_(0),
    vidNestedWPBitmap_(0),
    electronVeto_(0),
    mvaID_(0),
    mvaID_(0),
    pixelSeed_(0),
    genPartIdx_(0),
    genPartFlav_(0),
    cleanmask_(0)
  {}
  Float_t eCorr() const {return eCorr_;}
  Float_t energyErr() const {return energyErr_;}
  Float_t hoe() const {return hoe_;}
  Float_t mvaID() const {return mvaID_;}
  Float_t pfRelIso03() const {return pfRelIso03_;}
  Float_t pfRelIso03() const {return pfRelIso03_;}
  Float_t r9() const {return r9_;}
  Float_t sieie() const {return sieie_;}
  Int_t charge() const {return charge_;}
  Int_t cutBased() const {return cutBased_;}
  Int_t electronIdx() const {return electronIdx_;}
  Int_t jetIdx() const {return jetIdx_;}
  Int_t pdgId() const {return pdgId_;}
  Int_t vidNestedWPBitmap() const {return vidNestedWPBitmap_;}
  Bool_t electronVeto() const {return electronVeto_;}
  Bool_t mvaID() const {return mvaID_;}
  Bool_t mvaID() const {return mvaID_;}
  Bool_t pixelSeed() const {return pixelSeed_;}
  Int_t genPartIdx() const {return genPartIdx_;}
  UChar_t genPartFlav() const {return genPartFlav_;}
  UChar_t cleanmask() const {return cleanmask_;}
private:
  Float_t eCorr_;
  Float_t energyErr_;
  Float_t hoe_;
  Float_t mvaID_;
  Float_t pfRelIso03_;
  Float_t pfRelIso03_;
  Float_t r9_;
  Float_t sieie_;
  Int_t charge_;
  Int_t cutBased_;
  Int_t electronIdx_;
  Int_t jetIdx_;
  Int_t pdgId_;
  Int_t vidNestedWPBitmap_;
  Bool_t electronVeto_;
  Bool_t mvaID_;
  Bool_t mvaID_;
  Bool_t pixelSeed_;
  Int_t genPartIdx_;
  UChar_t genPartFlav_;
  UChar_t cleanmask_;
  void seteCorr(const Float_t value) {eCorr_ = value;}
  void setenergyErr(const Float_t value) {energyErr_ = value;}
  void sethoe(const Float_t value) {hoe_ = value;}
  void setmvaID(const Float_t value) {mvaID_ = value;}
  void setpfRelIso03(const Float_t value) {pfRelIso03_ = value;}
  void setpfRelIso03(const Float_t value) {pfRelIso03_ = value;}
  void setr9(const Float_t value) {r9_ = value;}
  void setsieie(const Float_t value) {sieie_ = value;}
  void setcharge(const Int_t value) {charge_ = value;}
  void setcutBased(const Int_t value) {cutBased_ = value;}
  void setelectronIdx(const Int_t value) {electronIdx_ = value;}
  void setjetIdx(const Int_t value) {jetIdx_ = value;}
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setvidNestedWPBitmap(const Int_t value) {vidNestedWPBitmap_ = value;}
  void setelectronVeto(const Bool_t value) {electronVeto_ = value;}
  void setmvaID(const Bool_t value) {mvaID_ = value;}
  void setmvaID(const Bool_t value) {mvaID_ = value;}
  void setpixelSeed(const Bool_t value) {pixelSeed_ = value;}
  void setgenPartIdx(const Int_t value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_t value) {genPartFlav_ = value;}
  void setcleanmask(const UChar_t value) {cleanmask_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Gendressedlepton: public TLorentzVector{
friend class URStreamer;
public:
//  Gendressedlepton(const Int_t &i_pdgId_):
//    
//  {}
  Gendressedlepton():
    TLorentzVector(),
    pdgId_(0)
  {}
  Int_t pdgId() const {return pdgId_;}
private:
  Int_t pdgId_;
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Softactivityjet: public TLorentzVector{
friend class URStreamer;
public:
//  Softactivityjet(const Float_t &i_SoftActivityJetHT_,const Float_t &i_SoftActivityJetHT10_,const Float_t &i_SoftActivityJetHT2_,const Float_t &i_SoftActivityJetHT5_,const Int_t &i_SoftActivityJetNjets10_,const Int_t &i_SoftActivityJetNjets2_,const Int_t &i_SoftActivityJetNjets5_):
//    
//  {}
  Softactivityjet():
    TLorentzVector(),
    SoftActivityJetHT_(0),
    SoftActivityJetHT10_(0),
    SoftActivityJetHT2_(0),
    SoftActivityJetHT5_(0),
    SoftActivityJetNjets10_(0),
    SoftActivityJetNjets2_(0),
    SoftActivityJetNjets5_(0)
  {}
  Float_t SoftActivityJetHT() const {return SoftActivityJetHT_;}
  Float_t SoftActivityJetHT10() const {return SoftActivityJetHT10_;}
  Float_t SoftActivityJetHT2() const {return SoftActivityJetHT2_;}
  Float_t SoftActivityJetHT5() const {return SoftActivityJetHT5_;}
  Int_t SoftActivityJetNjets10() const {return SoftActivityJetNjets10_;}
  Int_t SoftActivityJetNjets2() const {return SoftActivityJetNjets2_;}
  Int_t SoftActivityJetNjets5() const {return SoftActivityJetNjets5_;}
private:
  Float_t SoftActivityJetHT_;
  Float_t SoftActivityJetHT10_;
  Float_t SoftActivityJetHT2_;
  Float_t SoftActivityJetHT5_;
  Int_t SoftActivityJetNjets10_;
  Int_t SoftActivityJetNjets2_;
  Int_t SoftActivityJetNjets5_;
  void setSoftActivityJetHT(const Float_t value) {SoftActivityJetHT_ = value;}
  void setSoftActivityJetHT10(const Float_t value) {SoftActivityJetHT10_ = value;}
  void setSoftActivityJetHT2(const Float_t value) {SoftActivityJetHT2_ = value;}
  void setSoftActivityJetHT5(const Float_t value) {SoftActivityJetHT5_ = value;}
  void setSoftActivityJetNjets10(const Int_t value) {SoftActivityJetNjets10_ = value;}
  void setSoftActivityJetNjets2(const Int_t value) {SoftActivityJetNjets2_ = value;}
  void setSoftActivityJetNjets5(const Int_t value) {SoftActivityJetNjets5_ = value;}
  void setLotentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Subjet: public TLorentzVector{
friend class URStreamer;
public:
//  Subjet(const Float_t &i_btagCMVA_,const Float_t &i_btagCSVV2_,const Float_t &i_btagDeepB_,const Float_t &i_n2b1_,const Float_t &i_n3b1_,const Float_t &i_tau1_,const Float_t &i_tau2_,const Float_t &i_tau3_,const Float_t &i_tau4_):
//    
//  {}
  Subjet():
    TLorentzVector(),
    btagCMVA_(0),
    btagCSVV2_(0),
    btagDeepB_(0),
    n2b1_(0),
    n3b1_(0),
    tau1_(0),
    tau2_(0),
    tau3_(0),
    tau4_(0)
  {}
  Float_t btagCMVA() const {return btagCMVA_;}
  Float_t btagCSVV2() const {return btagCSVV2_;}
  Float_t btagDeepB() const {return btagDeepB_;}
  Float_t n2b1() const {return n2b1_;}
  Float_t n3b1() const {return n3b1_;}
  Float_t tau1() const {return tau1_;}
  Float_t tau2() const {return tau2_;}
  Float_t tau3() const {return tau3_;}
  Float_t tau4() const {return tau4_;}
private:
  Float_t btagCMVA_;
  Float_t btagCSVV2_;
  Float_t btagDeepB_;
  Float_t n2b1_;
  Float_t n3b1_;
  Float_t tau1_;
  Float_t tau2_;
  Float_t tau3_;
  Float_t tau4_;
  void setbtagCMVA(const Float_t value) {btagCMVA_ = value;}
  void setbtagCSVV2(const Float_t value) {btagCSVV2_ = value;}
  void setbtagDeepB(const Float_t value) {btagDeepB_ = value;}
  void setn2b1(const Float_t value) {n2b1_ = value;}
  void setn3b1(const Float_t value) {n3b1_ = value;}
  void settau1(const Float_t value) {tau1_ = value;}
  void settau2(const Float_t value) {tau2_ = value;}
  void settau3(const Float_t value) {tau3_ = value;}
  void settau4(const Float_t value) {tau4_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Tau: public TLorentzVector{
friend class URStreamer;
public:
//  Tau(const Float_t &i_chargedIso_,const Float_t &i_dxy_,const Float_t &i_dz_,const Float_t &i_footprintCorr_,const Float_t &i_leadTkDeltaEta_,const Float_t &i_leadTkDeltaPhi_,const Float_t &i_leadTkPtOverTauPt_,const Float_t &i_neutralIso_,const Float_t &i_photonsOutsideSignalCone_,const Float_t &i_puCorr_,const Float_t &i_rawAntiEle_,const Float_t &i_rawIso_,const Float_t &i_rawMVAnewDM_,const Float_t &i_rawMVAoldDM_,const Float_t &i_rawMVAoldDMdR03_,const Int_t &i_charge_,const Int_t &i_decayMode_,const Int_t &i_jetIdx_,const Int_t &i_rawAntiEleCat_,const UChar_t &i_idAntiEle_,const UChar_t &i_idAntiMu_,const Bool_t &i_idDecayMode_,const Bool_t &i_idDecayModeNewDMs_,const UChar_t &i_idMVAnewDM_,const UChar_t &i_idMVAoldDM_,const UChar_t &i_idMVAoldDMdR03_,const UChar_t &i_cleanmask_,const Int_t &i_genPartIdx_,const UChar_t &i_genPartFlav_):
//    
//  {}
  Tau():
    TLorentzVector(),
    chargedIso_(0),
    dxy_(0),
    dz_(0),
    footprintCorr_(0),
    leadTkDeltaEta_(0),
    leadTkDeltaPhi_(0),
    leadTkPtOverTauPt_(0),
    neutralIso_(0),
    photonsOutsideSignalCone_(0),
    puCorr_(0),
    rawAntiEle_(0),
    rawIso_(0),
    rawMVAnewDM_(0),
    rawMVAoldDM_(0),
    rawMVAoldDMdR03_(0),
    charge_(0),
    decayMode_(0),
    jetIdx_(0),
    rawAntiEleCat_(0),
    idAntiEle_(0),
    idAntiMu_(0),
    idDecayMode_(0),
    idDecayModeNewDMs_(0),
    idMVAnewDM_(0),
    idMVAoldDM_(0),
    idMVAoldDMdR03_(0),
    cleanmask_(0),
    genPartIdx_(0),
    genPartFlav_(0)
  {}
  Float_t chargedIso() const {return chargedIso_;}
  Float_t dxy() const {return dxy_;}
  Float_t dz() const {return dz_;}
  Float_t footprintCorr() const {return footprintCorr_;}
  Float_t leadTkDeltaEta() const {return leadTkDeltaEta_;}
  Float_t leadTkDeltaPhi() const {return leadTkDeltaPhi_;}
  Float_t leadTkPtOverTauPt() const {return leadTkPtOverTauPt_;}
  Float_t neutralIso() const {return neutralIso_;}
  Float_t photonsOutsideSignalCone() const {return photonsOutsideSignalCone_;}
  Float_t puCorr() const {return puCorr_;}
  Float_t rawAntiEle() const {return rawAntiEle_;}
  Float_t rawIso() const {return rawIso_;}
  Float_t rawMVAnewDM() const {return rawMVAnewDM_;}
  Float_t rawMVAoldDM() const {return rawMVAoldDM_;}
  Float_t rawMVAoldDMdR03() const {return rawMVAoldDMdR03_;}
  Int_t charge() const {return charge_;}
  Int_t decayMode() const {return decayMode_;}
  Int_t jetIdx() const {return jetIdx_;}
  Int_t rawAntiEleCat() const {return rawAntiEleCat_;}
  UChar_t idAntiEle() const {return idAntiEle_;}
  UChar_t idAntiMu() const {return idAntiMu_;}
  Bool_t idDecayMode() const {return idDecayMode_;}
  Bool_t idDecayModeNewDMs() const {return idDecayModeNewDMs_;}
  UChar_t idMVAnewDM() const {return idMVAnewDM_;}
  UChar_t idMVAoldDM() const {return idMVAoldDM_;}
  UChar_t idMVAoldDMdR03() const {return idMVAoldDMdR03_;}
  UChar_t cleanmask() const {return cleanmask_;}
  Int_t genPartIdx() const {return genPartIdx_;}
  UChar_t genPartFlav() const {return genPartFlav_;}
private:
  Float_t chargedIso_;
  Float_t dxy_;
  Float_t dz_;
  Float_t footprintCorr_;
  Float_t leadTkDeltaEta_;
  Float_t leadTkDeltaPhi_;
  Float_t leadTkPtOverTauPt_;
  Float_t neutralIso_;
  Float_t photonsOutsideSignalCone_;
  Float_t puCorr_;
  Float_t rawAntiEle_;
  Float_t rawIso_;
  Float_t rawMVAnewDM_;
  Float_t rawMVAoldDM_;
  Float_t rawMVAoldDMdR03_;
  Int_t charge_;
  Int_t decayMode_;
  Int_t jetIdx_;
  Int_t rawAntiEleCat_;
  UChar_t idAntiEle_;
  UChar_t idAntiMu_;
  Bool_t idDecayMode_;
  Bool_t idDecayModeNewDMs_;
  UChar_t idMVAnewDM_;
  UChar_t idMVAoldDM_;
  UChar_t idMVAoldDMdR03_;
  UChar_t cleanmask_;
  Int_t genPartIdx_;
  UChar_t genPartFlav_;
  void setchargedIso(const Float_t value) {chargedIso_ = value;}
  void setdxy(const Float_t value) {dxy_ = value;}
  void setdz(const Float_t value) {dz_ = value;}
  void setfootprintCorr(const Float_t value) {footprintCorr_ = value;}
  void setleadTkDeltaEta(const Float_t value) {leadTkDeltaEta_ = value;}
  void setleadTkDeltaPhi(const Float_t value) {leadTkDeltaPhi_ = value;}
  void setleadTkPtOverTauPt(const Float_t value) {leadTkPtOverTauPt_ = value;}
  void setneutralIso(const Float_t value) {neutralIso_ = value;}
  void setphotonsOutsideSignalCone(const Float_t value) {photonsOutsideSignalCone_ = value;}
  void setpuCorr(const Float_t value) {puCorr_ = value;}
  void setrawAntiEle(const Float_t value) {rawAntiEle_ = value;}
  void setrawIso(const Float_t value) {rawIso_ = value;}
  void setrawMVAnewDM(const Float_t value) {rawMVAnewDM_ = value;}
  void setrawMVAoldDM(const Float_t value) {rawMVAoldDM_ = value;}
  void setrawMVAoldDMdR03(const Float_t value) {rawMVAoldDMdR03_ = value;}
  void setcharge(const Int_t value) {charge_ = value;}
  void setdecayMode(const Int_t value) {decayMode_ = value;}
  void setjetIdx(const Int_t value) {jetIdx_ = value;}
  void setrawAntiEleCat(const Int_t value) {rawAntiEleCat_ = value;}
  void setidAntiEle(const UChar_t value) {idAntiEle_ = value;}
  void setidAntiMu(const UChar_t value) {idAntiMu_ = value;}
  void setidDecayMode(const Bool_t value) {idDecayMode_ = value;}
  void setidDecayModeNewDMs(const Bool_t value) {idDecayModeNewDMs_ = value;}
  void setidMVAnewDM(const UChar_t value) {idMVAnewDM_ = value;}
  void setidMVAoldDM(const UChar_t value) {idMVAoldDM_ = value;}
  void setidMVAoldDMdR03(const UChar_t value) {idMVAoldDMdR03_ = value;}
  void setcleanmask(const UChar_t value) {cleanmask_ = value;}
  void setgenPartIdx(const Int_t value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_t value) {genPartFlav_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Trigobj: public TLorentzVector{
friend class URStreamer;
public:
//  Trigobj(const Float_t &i_l1pt_,const Float_t &i_l1pt_,const Float_t &i_l2pt_,const Int_t &i_id_,const Int_t &i_l1iso_,const Int_t &i_l1charge_,const Int_t &i_filterBits_):
//    
//  {}
  Trigobj():
    TLorentzVector(),
    l1pt_(0),
    l1pt_(0),
    l2pt_(0),
    id_(0),
    l1iso_(0),
    l1charge_(0),
    filterBits_(0)
  {}
  Float_t l1pt() const {return l1pt_;}
  Float_t l1pt() const {return l1pt_;}
  Float_t l2pt() const {return l2pt_;}
  Int_t id() const {return id_;}
  Int_t l1iso() const {return l1iso_;}
  Int_t l1charge() const {return l1charge_;}
  Int_t filterBits() const {return filterBits_;}
private:
  Float_t l1pt_;
  Float_t l1pt_;
  Float_t l2pt_;
  Int_t id_;
  Int_t l1iso_;
  Int_t l1charge_;
  Int_t filterBits_;
  void setl1pt(const Float_t value) {l1pt_ = value;}
  void setl1pt(const Float_t value) {l1pt_ = value;}
  void setl2pt(const Float_t value) {l2pt_ = value;}
  void setid(const Int_t value) {id_ = value;}
  void setl1iso(const Int_t value) {l1iso_ = value;}
  void setl1charge(const Int_t value) {l1charge_ = value;}
  void setfilterBits(const Int_t value) {filterBits_ = value;}
  void setLotentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Otherpv{
friend class URStreamer;
public:
//  Otherpv(const Float_t &i_z_):
//    
//  {}
  Otherpv():
    z_(0)
  {}
  Float_t z() const {return z_;}
private:
  Float_t z_;
  void setz(const Float_t value) {z_ = value;}
};

class Sv: public TLorentzVector{
friend class URStreamer;
public:
//  Sv(const Float_t &i_dlen_,const Float_t &i_dlenSig_,const Float_t &i_pAngle_,const Float_t &i_chi2_,const Float_t &i_ndof_,const Float_t &i_x_,const Float_t &i_y_,const Float_t &i_z_):
//    
//  {}
  Sv():
    TLorentzVector(),
    dlen_(0),
    dlenSig_(0),
    pAngle_(0),
    chi2_(0),
    ndof_(0),
    x_(0),
    y_(0),
    z_(0)
  {}
  Float_t dlen() const {return dlen_;}
  Float_t dlenSig() const {return dlenSig_;}
  Float_t pAngle() const {return pAngle_;}
  Float_t chi2() const {return chi2_;}
  Float_t ndof() const {return ndof_;}
  Float_t x() const {return x_;}
  Float_t y() const {return y_;}
  Float_t z() const {return z_;}
private:
  Float_t dlen_;
  Float_t dlenSig_;
  Float_t pAngle_;
  Float_t chi2_;
  Float_t ndof_;
  Float_t x_;
  Float_t y_;
  Float_t z_;
  void setdlen(const Float_t value) {dlen_ = value;}
  void setdlenSig(const Float_t value) {dlenSig_ = value;}
  void setpAngle(const Float_t value) {pAngle_ = value;}
  void setchi2(const Float_t value) {chi2_ = value;}
  void setndof(const Float_t value) {ndof_ = value;}
  void setx(const Float_t value) {x_ = value;}
  void sety(const Float_t value) {y_ = value;}
  void setz(const Float_t value) {z_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};


class URStreamer{
public:
  UInt_t run;
  UInt_t luminosityBlock;
  ULong64_t event;
  UInt_t nElectron;
  UInt_t nFatJet;
  UInt_t nGenJetAK8;
  UInt_t nGenJet;
  UInt_t nGenPart;
  UInt_t nGenVisTau;
  Float_t genWeight;
  UInt_t nLHEPdfWeight;
  Float_t LHEPdfWeight;
  UInt_t nLHEScaleWeight;
  Float_t LHEScaleWeight;
  UInt_t nJet;
  UInt_t nMuon;
  UInt_t nPhoton;
  Float_t fixedGridRhoFastjetAll;
  Float_t fixedGridRhoFastjetCentralCalo;
  Float_t fixedGridRhoFastjetCentralNeutral;
  UInt_t nGenDressedLepton;
  UInt_t nSoftActivityJet;
  Float_t SoftActivityJetHT;
  Float_t SoftActivityJetHT10;
  Float_t SoftActivityJetHT2;
  Float_t SoftActivityJetHT5;
  Int_t SoftActivityJetNjets10;
  Int_t SoftActivityJetNjets2;
  Int_t SoftActivityJetNjets5;
  UInt_t nSubJet;
  UInt_t nTau;
  UInt_t nTrigObj;
  Int_t genTtbarId;
  UInt_t nOtherPV;
  UInt_t nSV;
  Bool_t HLTriggerFirstPath;
  Bool_t HLTriggerFinalPath;

  URStreamer(TTree *tree):
    run(0),
    luminosityBlock(0),
    event(0),
    CaloMET_phi_(0),
    CaloMET_pt_(0),
    CaloMET_sumEt_(0),
    nElectron(0),
    Electron_deltaEtaSC_(0),
    Electron_dr03EcalRecHitSumEt_(0),
    Electron_dr03HcalDepth1TowerSumEt_(0),
    Electron_dr03TkSumPt_(0),
    Electron_dxy_(0),
    Electron_dxyErr_(0),
    Electron_dz_(0),
    Electron_dzErr_(0),
    Electron_eCorr_(0),
    Electron_eInvMinusPInv_(0),
    Electron_energyErr_(0),
    Electron_eta_(0),
    Electron_hoe_(0),
    Electron_ip3d_(0),
    Electron_mass_(0),
    Electron_miniPFRelIso_all_(0),
    Electron_miniPFRelIso_chg_(0),
    Electron_mvaSpring16GP_(0),
    Electron_mvaSpring16HZZ_(0),
    Electron_pfRelIso03_all_(0),
    Electron_pfRelIso03_chg_(0),
    Electron_phi_(0),
    Electron_pt_(0),
    Electron_r9_(0),
    Electron_sieie_(0),
    Electron_sip3d_(0),
    Electron_mvaTTH_(0),
    Electron_charge_(0),
    Electron_cutBased_(0),
    Electron_cutBased_HLTPreSel_(0),
    Electron_jetIdx_(0),
    Electron_pdgId_(0),
    Electron_photonIdx_(0),
    Electron_tightCharge_(0),
    Electron_vidNestedWPBitmap_(0),
    Electron_convVeto_(0),
    Electron_cutBased_HEEP_(0),
    Electron_isPFcand_(0),
    Electron_lostHits_(0),
    Electron_mvaSpring16GP_WP80_(0),
    Electron_mvaSpring16GP_WP90_(0),
    Electron_mvaSpring16HZZ_WPL_(0),
    nFatJet(0),
    FatJet_area_(0),
    FatJet_btagCMVA_(0),
    FatJet_btagCSVV2_(0),
    FatJet_btagDeepB_(0),
    FatJet_btagHbb_(0),
    FatJet_eta_(0),
    FatJet_mass_(0),
    FatJet_msoftdrop_(0),
    FatJet_msoftdrop_chs_(0),
    FatJet_n2b1_(0),
    FatJet_n3b1_(0),
    FatJet_phi_(0),
    FatJet_pt_(0),
    FatJet_tau1_(0),
    FatJet_tau2_(0),
    FatJet_tau3_(0),
    FatJet_tau4_(0),
    FatJet_jetId_(0),
    FatJet_subJetIdx1_(0),
    FatJet_subJetIdx2_(0),
    nGenJetAK8(0),
    GenJetAK8_eta_(0),
    GenJetAK8_mass_(0),
    GenJetAK8_phi_(0),
    GenJetAK8_pt_(0),
    nGenJet(0),
    GenJet_eta_(0),
    GenJet_mass_(0),
    GenJet_phi_(0),
    GenJet_pt_(0),
    nGenPart(0),
    GenPart_eta_(0),
    GenPart_mass_(0),
    GenPart_phi_(0),
    GenPart_pt_(0),
    GenPart_genPartIdxMother_(0),
    GenPart_pdgId_(0),
    GenPart_status_(0),
    GenPart_statusFlags_(0),
    Generator_binvar_(0),
    Generator_scalePDF_(0),
    Generator_weight_(0),
    Generator_x1_(0),
    Generator_x2_(0),
    Generator_xpdf1_(0),
    Generator_xpdf2_(0),
    Generator_id1_(0),
    Generator_id2_(0),
    nGenVisTau(0),
    GenVisTau_eta_(0),
    GenVisTau_mass_(0),
    GenVisTau_phi_(0),
    GenVisTau_pt_(0),
    GenVisTau_charge_(0),
    GenVisTau_genPartIdxMother_(0),
    GenVisTau_status_(0),
    genWeight(0),
    LHEWeight_originalXWGTUP_(0),
    nLHEPdfWeight(0),
    LHEPdfWeight(0),
    nLHEScaleWeight(0),
    LHEScaleWeight(0),
    nJet(0),
    Jet_area_(0),
    Jet_btagCMVA_(0),
    Jet_btagCSVV2_(0),
    Jet_btagDeepB_(0),
    Jet_btagDeepC_(0),
    Jet_chEmEF_(0),
    Jet_chHEF_(0),
    Jet_eta_(0),
    Jet_mass_(0),
    Jet_neEmEF_(0),
    Jet_neHEF_(0),
    Jet_phi_(0),
    Jet_pt_(0),
    Jet_qgl_(0),
    Jet_rawFactor_(0),
    Jet_bReg_(0),
    Jet_electronIdx1_(0),
    Jet_electronIdx2_(0),
    Jet_jetId_(0),
    Jet_muonIdx1_(0),
    Jet_muonIdx2_(0),
    Jet_nConstituents_(0),
    Jet_nElectrons_(0),
    Jet_nMuons_(0),
    Jet_puId_(0),
    LHE_HT_(0),
    LHE_HTIncoming_(0),
    LHE_Vpt_(0),
    LHE_Njets_(0),
    LHE_Nb_(0),
    LHE_Nc_(0),
    LHE_Nuds_(0),
    LHE_Nglu_(0),
    LHE_NpNLO_(0),
    LHE_NpLO_(0),
    GenMET_phi_(0),
    GenMET_pt_(0),
    MET_MetUnclustEnUpDeltaX_(0),
    MET_MetUnclustEnUpDeltaY_(0),
    MET_covXX_(0),
    MET_covXY_(0),
    MET_covYY_(0),
    MET_phi_(0),
    MET_pt_(0),
    MET_significance_(0),
    MET_sumEt_(0),
    nMuon(0),
    Muon_dxy_(0),
    Muon_dxyErr_(0),
    Muon_dz_(0),
    Muon_dzErr_(0),
    Muon_eta_(0),
    Muon_ip3d_(0),
    Muon_mass_(0),
    Muon_miniPFRelIso_all_(0),
    Muon_miniPFRelIso_chg_(0),
    Muon_pfRelIso03_all_(0),
    Muon_pfRelIso03_chg_(0),
    Muon_pfRelIso04_all_(0),
    Muon_phi_(0),
    Muon_pt_(0),
    Muon_ptErr_(0),
    Muon_segmentComp_(0),
    Muon_sip3d_(0),
    Muon_mvaTTH_(0),
    Muon_charge_(0),
    Muon_jetIdx_(0),
    Muon_nStations_(0),
    Muon_nTrackerLayers_(0),
    Muon_pdgId_(0),
    Muon_tightCharge_(0),
    Muon_highPtId_(0),
    Muon_isPFcand_(0),
    Muon_mediumId_(0),
    Muon_softId_(0),
    Muon_tightId_(0),
    nPhoton(0),
    Photon_eCorr_(0),
    Photon_energyErr_(0),
    Photon_eta_(0),
    Photon_hoe_(0),
    Photon_mass_(0),
    Photon_mvaID_(0),
    Photon_pfRelIso03_all_(0),
    Photon_pfRelIso03_chg_(0),
    Photon_phi_(0),
    Photon_pt_(0),
    Photon_r9_(0),
    Photon_sieie_(0),
    Photon_charge_(0),
    Photon_cutBased_(0),
    Photon_electronIdx_(0),
    Photon_jetIdx_(0),
    Photon_pdgId_(0),
    Photon_vidNestedWPBitmap_(0),
    Photon_electronVeto_(0),
    Photon_mvaID_WP80_(0),
    Photon_mvaID_WP90_(0),
    Photon_pixelSeed_(0),
    Pileup_nTrueInt_(0),
    Pileup_nPU_(0),
    Pileup_sumEOOT_(0),
    Pileup_sumLOOT_(0),
    PuppiMET_phi_(0),
    PuppiMET_pt_(0),
    PuppiMET_sumEt_(0),
    RawMET_phi_(0),
    RawMET_pt_(0),
    RawMET_sumEt_(0),
    fixedGridRhoFastjetAll(0),
    fixedGridRhoFastjetCentralCalo(0),
    fixedGridRhoFastjetCentralNeutral(0),
    nGenDressedLepton(0),
    GenDressedLepton_eta_(0),
    GenDressedLepton_mass_(0),
    GenDressedLepton_phi_(0),
    GenDressedLepton_pt_(0),
    GenDressedLepton_pdgId_(0),
    nSoftActivityJet(0),
    SoftActivityJet_eta_(0),
    SoftActivityJet_phi_(0),
    SoftActivityJet_pt_(0),
    SoftActivityJetHT(0),
    SoftActivityJetHT10(0),
    SoftActivityJetHT2(0),
    SoftActivityJetHT5(0),
    SoftActivityJetNjets10(0),
    SoftActivityJetNjets2(0),
    SoftActivityJetNjets5(0),
    nSubJet(0),
    SubJet_btagCMVA_(0),
    SubJet_btagCSVV2_(0),
    SubJet_btagDeepB_(0),
    SubJet_eta_(0),
    SubJet_mass_(0),
    SubJet_n2b1_(0),
    SubJet_n3b1_(0),
    SubJet_phi_(0),
    SubJet_pt_(0),
    SubJet_tau1_(0),
    SubJet_tau2_(0),
    SubJet_tau3_(0),
    SubJet_tau4_(0),
    nTau(0),
    Tau_chargedIso_(0),
    Tau_dxy_(0),
    Tau_dz_(0),
    Tau_eta_(0),
    Tau_footprintCorr_(0),
    Tau_leadTkDeltaEta_(0),
    Tau_leadTkDeltaPhi_(0),
    Tau_leadTkPtOverTauPt_(0),
    Tau_mass_(0),
    Tau_neutralIso_(0),
    Tau_phi_(0),
    Tau_photonsOutsideSignalCone_(0),
    Tau_pt_(0),
    Tau_puCorr_(0),
    Tau_rawAntiEle_(0),
    Tau_rawIso_(0),
    Tau_rawMVAnewDM_(0),
    Tau_rawMVAoldDM_(0),
    Tau_rawMVAoldDMdR03_(0),
    Tau_charge_(0),
    Tau_decayMode_(0),
    Tau_jetIdx_(0),
    Tau_rawAntiEleCat_(0),
    Tau_idAntiEle_(0),
    Tau_idAntiMu_(0),
    Tau_idDecayMode_(0),
    Tau_idDecayModeNewDMs_(0),
    Tau_idMVAnewDM_(0),
    Tau_idMVAoldDM_(0),
    Tau_idMVAoldDMdR03_(0),
    TkMET_phi_(0),
    TkMET_pt_(0),
    TkMET_sumEt_(0),
    nTrigObj(0),
    TrigObj_pt_(0),
    TrigObj_eta_(0),
    TrigObj_phi_(0),
    TrigObj_l1pt_(0),
    TrigObj_l1pt_2_(0),
    TrigObj_l2pt_(0),
    TrigObj_id_(0),
    TrigObj_l1iso_(0),
    TrigObj_l1charge_(0),
    TrigObj_filterBits_(0),
    genTtbarId(0),
    nOtherPV(0),
    OtherPV_z_(0),
    PV_ndof_(0),
    PV_x_(0),
    PV_y_(0),
    PV_z_(0),
    PV_chi2_(0),
    PV_score_(0),
    PV_npvs_(0),
    PV_npvsGood_(0),
    nSV(0),
    SV_dlen_(0),
    SV_dlenSig_(0),
    SV_pAngle_(0),
    Electron_genPartIdx_(0),
    Electron_genPartFlav_(0),
    GenJetAK8_partonFlavour_(0),
    GenJetAK8_hadronFlavour_(0),
    GenJet_partonFlavour_(0),
    GenJet_hadronFlavour_(0),
    Jet_genJetIdx_(0),
    Jet_hadronFlavour_(0),
    Jet_partonFlavour_(0),
    Muon_genPartIdx_(0),
    Muon_genPartFlav_(0),
    Photon_genPartIdx_(0),
    Photon_genPartFlav_(0),
    MET_fiducialGenPhi_(0),
    MET_fiducialGenPt_(0),
    Electron_cleanmask_(0),
    Jet_cleanmask_(0),
    Muon_cleanmask_(0),
    Photon_cleanmask_(0),
    Tau_cleanmask_(0),
    SV_chi2_(0),
    SV_eta_(0),
    SV_mass_(0),
    SV_ndof_(0),
    SV_phi_(0),
    SV_pt_(0),
    SV_x_(0),
    SV_y_(0),
    SV_z_(0),
    Tau_genPartIdx_(0),
    Tau_genPartFlav_(0),
    L1simulation_step_(0),
    HLTriggerFirstPath(0),
    HLT_AK8PFJet360_TrimMass30_(0),
    HLT_AK8PFJet400_TrimMass30_(0),
    HLT_AK8PFHT750_TrimMass50_(0),
    HLT_AK8PFHT800_TrimMass50_(0),
    HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_(0),
    HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_(0),
    HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_(0),
    HLT_AK8DiPFJet300_200_TrimMass30_(0),
    HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_(0),
    HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_(0),
    HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_(0),
    HLT_AK8DiPFJet280_200_TrimMass30_(0),
    HLT_AK8DiPFJet250_200_TrimMass30_(0),
    HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_(0),
    HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_(0),
    HLT_CaloJet260_(0),
    HLT_CaloJet500_NoJetID_(0),
    HLT_Dimuon13_PsiPrime_(0),
    HLT_Dimuon13_Upsilon_(0),
    HLT_Dimuon20_Jpsi_(0),
    HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_(0),
    HLT_DoubleEle25_CaloIdL_GsfTrkIdVL_(0),
    HLT_DoubleEle33_CaloIdL_(0),
    HLT_DoubleEle33_CaloIdL_MW_(0),
    HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_(0),
    HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_(0),
    HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_(0),
    HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_(0),
    HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_(0),
    HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_(0),
    HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_(0),
    HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_(0),
    HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_(0),
    HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_(0),
    HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_(0),
    HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_(0),
    HLT_DoubleMu33NoFiltersNoVtx_(0),
    HLT_DoubleMu38NoFiltersNoVtx_(0),
    HLT_DoubleMu23NoFiltersNoVtxDisplaced_(0),
    HLT_DoubleMu28NoFiltersNoVtxDisplaced_(0),
    HLT_DoubleMu0_(0),
    HLT_DoubleMu4_3_Bs_(0),
    HLT_DoubleMu4_3_Jpsi_Displaced_(0),
    HLT_DoubleMu4_JpsiTrk_Displaced_(0),
    HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_(0),
    HLT_DoubleMu3_Trk_Tau3mu_(0),
    HLT_DoubleMu4_PsiPrimeTrk_Displaced_(0),
    HLT_Mu7p5_L2Mu2_Jpsi_(0),
    HLT_Mu7p5_L2Mu2_Upsilon_(0),
    HLT_Mu7p5_Track2_Jpsi_(0),
    HLT_Mu7p5_Track3p5_Jpsi_(0),
    HLT_Mu7p5_Track7_Jpsi_(0),
    HLT_Mu7p5_Track2_Upsilon_(0),
    HLT_Mu7p5_Track3p5_Upsilon_(0),
    HLT_Mu7p5_Track7_Upsilon_(0),
    HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_(0),
    HLT_Dimuon0er16_Jpsi_NoVertexing_(0),
    HLT_Dimuon6_Jpsi_NoVertexing_(0),
    HLT_Photon150_(0),
    HLT_Photon90_CaloIdL_HT300_(0),
    HLT_HT250_CaloMET70_(0),
    HLT_DoublePhoton60_(0),
    HLT_DoublePhoton85_(0),
    HLT_Ele17_Ele8_Gsf_(0),
    HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28_(0),
    HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29_(0),
    HLT_Ele22_eta2p1_WPLoose_Gsf_(0),
    HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_(0),
    HLT_Ele23_WPLoose_Gsf_(0),
    HLT_Ele23_WPLoose_Gsf_WHbbBoost_(0),
    HLT_Ele24_eta2p1_WPLoose_Gsf_(0),
    HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_(0),
    HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_(0),
    HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_(0),
    HLT_Ele25_WPTight_Gsf_(0),
    HLT_Ele25_eta2p1_WPLoose_Gsf_(0),
    HLT_Ele25_eta2p1_WPTight_Gsf_(0),
    HLT_Ele27_WPLoose_Gsf_(0),
    HLT_Ele27_WPLoose_Gsf_WHbbBoost_(0),
    HLT_Ele27_WPTight_Gsf_(0),
    HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_(0),
    HLT_Ele27_eta2p1_WPLoose_Gsf_(0),
    HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_(0),
    HLT_Ele27_eta2p1_WPTight_Gsf_(0),
    HLT_Ele30_WPTight_Gsf_(0),
    HLT_Ele30_eta2p1_WPLoose_Gsf_(0),
    HLT_Ele30_eta2p1_WPTight_Gsf_(0),
    HLT_Ele32_WPTight_Gsf_(0),
    HLT_Ele32_eta2p1_WPLoose_Gsf_(0),
    HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_(0),
    HLT_Ele32_eta2p1_WPTight_Gsf_(0),
    HLT_Ele35_WPLoose_Gsf_(0),
    HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_(0),
    HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_(0),
    HLT_Ele45_WPLoose_Gsf_(0),
    HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded_(0),
    HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_(0),
    HLT_Ele105_CaloIdVT_GsfTrkIdT_(0),
    HLT_Ele30WP60_SC4_Mass55_(0),
    HLT_Ele30WP60_Ele8_Mass55_(0),
    HLT_HT200_(0),
    HLT_HT275_(0),
    HLT_HT325_(0),
    HLT_HT425_(0),
    HLT_HT575_(0),
    HLT_HT410to430_(0),
    HLT_HT430to450_(0),
    HLT_HT450to470_(0),
    HLT_HT470to500_(0),
    HLT_HT500to550_(0),
    HLT_HT550to650_(0),
    HLT_HT650_(0),
    HLT_Mu16_eta2p1_MET30_(0),
    HLT_IsoMu16_eta2p1_MET30_(0),
    HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1_(0),
    HLT_IsoMu17_eta2p1_(0),
    HLT_IsoMu17_eta2p1_LooseIsoPFTau20_(0),
    HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_(0),
    HLT_DoubleIsoMu17_eta2p1_(0),
    HLT_DoubleIsoMu17_eta2p1_noDzCut_(0),
    HLT_IsoMu18_(0),
    HLT_IsoMu19_eta2p1_LooseIsoPFTau20_(0),
    HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_(0),
    HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_(0),
    HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_(0),
    HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_(0),
    HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_(0),
    HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_(0),
    HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_(0),
    HLT_IsoMu20_(0),
    HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_(0),
    HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1_(0),
    HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_(0),
    HLT_IsoMu22_(0),
    HLT_IsoMu22_eta2p1_(0),
    HLT_IsoMu24_(0),
    HLT_IsoMu27_(0),
    HLT_IsoTkMu18_(0),
    HLT_IsoTkMu20_(0),
    HLT_IsoTkMu22_(0),
    HLT_IsoTkMu22_eta2p1_(0),
    HLT_IsoTkMu24_(0),
    HLT_IsoTkMu27_(0),
    HLT_JetE30_NoBPTX3BX_(0),
    HLT_JetE30_NoBPTX_(0),
    HLT_JetE50_NoBPTX3BX_(0),
    HLT_JetE70_NoBPTX3BX_(0),
    HLT_L1SingleMu18_(0),
    HLT_L2Mu10_(0),
    HLT_L1SingleMuOpen_(0),
    HLT_L1SingleMuOpen_DT_(0),
    HLT_L2DoubleMu23_NoVertex_(0),
    HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_(0),
    HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10_(0),
    HLT_L2Mu10_NoVertex_NoBPTX3BX_(0),
    HLT_L2Mu10_NoVertex_NoBPTX_(0),
    HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_(0),
    HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_(0),
    HLT_LooseIsoPFTau50_Trk30_eta2p1_(0),
    HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_(0),
    HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_(0),
    HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_(0),
    HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_(0),
    HLT_PFTau120_eta2p1_(0),
    HLT_PFTau140_eta2p1_(0),
    HLT_VLooseIsoPFTau120_Trk50_eta2p1_(0),
    HLT_VLooseIsoPFTau140_Trk50_eta2p1_(0),
    HLT_Mu17_Mu8_(0),
    HLT_Mu17_Mu8_DZ_(0),
    HLT_Mu17_Mu8_SameSign_(0),
    HLT_Mu17_Mu8_SameSign_DZ_(0),
    HLT_Mu20_Mu10_(0),
    HLT_Mu20_Mu10_DZ_(0),
    HLT_Mu20_Mu10_SameSign_(0),
    HLT_Mu20_Mu10_SameSign_DZ_(0),
    HLT_Mu17_TkMu8_DZ_(0),
    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_(0),
    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_(0),
    HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_(0),
    HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_(0),
    HLT_Mu25_TkMu0_dEta18_Onia_(0),
    HLT_Mu27_TkMu8_(0),
    HLT_Mu30_TkMu11_(0),
    HLT_Mu30_eta2p1_PFJet150_PFJet50_(0),
    HLT_Mu40_TkMu11_(0),
    HLT_Mu40_eta2p1_PFJet200_PFJet50_(0),
    HLT_Mu20_(0),
    HLT_TkMu17_(0),
    HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_(0),
    HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_(0),
    HLT_TkMu20_(0),
    HLT_Mu24_eta2p1_(0),
    HLT_TkMu24_eta2p1_(0),
    HLT_Mu27_(0),
    HLT_TkMu27_(0),
    HLT_Mu45_eta2p1_(0),
    HLT_Mu50_(0),
    HLT_TkMu50_(0),
    HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_(0),
    HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL_(0),
    HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL_(0),
    HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL_(0),
    HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL_(0),
    HLT_DoubleMu18NoFiltersNoVtx_(0),
    HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight_(0),
    HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose_(0),
    HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose_(0),
    HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight_(0),
    HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose_(0),
    HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose_(0),
    HLT_Mu28NoFiltersNoVtx_CentralCaloJet40_(0),
    HLT_PFHT300_PFMET100_(0),
    HLT_PFHT300_PFMET110_(0),
    HLT_PFHT550_4JetPt50_(0),
    HLT_PFHT650_4JetPt50_(0),
    HLT_PFHT750_4JetPt50_(0),
    HLT_PFHT750_4JetPt70_(0),
    HLT_PFHT750_4JetPt80_(0),
    HLT_PFHT800_4JetPt50_(0),
    HLT_PFHT850_4JetPt50_(0),
    HLT_PFJet15_NoCaloMatched_(0),
    HLT_PFJet25_NoCaloMatched_(0),
    HLT_DiPFJet15_NoCaloMatched_(0),
    HLT_DiPFJet25_NoCaloMatched_(0),
    HLT_DiPFJet15_FBEta3_NoCaloMatched_(0),
    HLT_DiPFJet25_FBEta3_NoCaloMatched_(0),
    HLT_DiPFJetAve15_HFJEC_(0),
    HLT_DiPFJetAve25_HFJEC_(0),
    HLT_DiPFJetAve35_HFJEC_(0),
    HLT_AK8PFJet40_(0),
    HLT_AK8PFJet60_(0),
    HLT_AK8PFJet80_(0),
    HLT_AK8PFJet140_(0),
    HLT_AK8PFJet200_(0),
    HLT_AK8PFJet260_(0),
    HLT_AK8PFJet320_(0),
    HLT_AK8PFJet400_(0),
    HLT_AK8PFJet450_(0),
    HLT_AK8PFJet500_(0),
    HLT_PFJet40_(0),
    HLT_PFJet60_(0),
    HLT_PFJet80_(0),
    HLT_PFJet140_(0),
    HLT_PFJet200_(0),
    HLT_PFJet260_(0),
    HLT_PFJet320_(0),
    HLT_PFJet400_(0),
    HLT_PFJet450_(0),
    HLT_PFJet500_(0),
    HLT_DiPFJetAve40_(0),
    HLT_DiPFJetAve60_(0),
    HLT_DiPFJetAve80_(0),
    HLT_DiPFJetAve140_(0),
    HLT_DiPFJetAve200_(0),
    HLT_DiPFJetAve260_(0),
    HLT_DiPFJetAve320_(0),
    HLT_DiPFJetAve400_(0),
    HLT_DiPFJetAve500_(0),
    HLT_DiPFJetAve60_HFJEC_(0),
    HLT_DiPFJetAve80_HFJEC_(0),
    HLT_DiPFJetAve100_HFJEC_(0),
    HLT_DiPFJetAve160_HFJEC_(0),
    HLT_DiPFJetAve220_HFJEC_(0),
    HLT_DiPFJetAve300_HFJEC_(0),
    HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140_(0),
    HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80_(0),
    HLT_DiCentralPFJet170_(0),
    HLT_SingleCentralPFJet170_CFMax0p1_(0),
    HLT_DiCentralPFJet170_CFMax0p1_(0),
    HLT_DiCentralPFJet220_CFMax0p3_(0),
    HLT_DiCentralPFJet330_CFMax0p5_(0),
    HLT_DiCentralPFJet430_(0),
    HLT_PFHT125_(0),
    HLT_PFHT200_(0),
    HLT_PFHT250_(0),
    HLT_PFHT300_(0),
    HLT_PFHT350_(0),
    HLT_PFHT400_(0),
    HLT_PFHT475_(0),
    HLT_PFHT600_(0),
    HLT_PFHT650_(0),
    HLT_PFHT800_(0),
    HLT_PFHT900_(0),
    HLT_PFHT200_PFAlphaT0p51_(0),
    HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57_(0),
    HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_(0),
    HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55_(0),
    HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_(0),
    HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53_(0),
    HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_(0),
    HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52_(0),
    HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53_(0),
    HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51_(0),
    HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52_(0),
    HLT_MET60_IsoTrk35_Loose_(0),
    HLT_MET75_IsoTrk50_(0),
    HLT_MET90_IsoTrk50_(0),
    HLT_PFMET120_BTagCSV_p067_(0),
    HLT_PFMET120_Mu5_(0),
    HLT_PFMET170_NotCleaned_(0),
    HLT_PFMET170_NoiseCleaned_(0),
    HLT_PFMET170_HBHECleaned_(0),
    HLT_PFMET170_JetIdCleaned_(0),
    HLT_PFMET170_BeamHaloCleaned_(0),
    HLT_PFMET170_HBHE_BeamHaloCleaned_(0),
    HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned_(0),
    HLT_PFMET90_PFMHT90_IDTight_(0),
    HLT_PFMET100_PFMHT100_IDTight_(0),
    HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned_(0),
    HLT_PFMET110_PFMHT110_IDTight_(0),
    HLT_PFMET120_PFMHT120_IDTight_(0),
    HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_(0),
    HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_(0),
    HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200_(0),
    HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460_(0),
    HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_(0),
    HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_(0),
    HLT_QuadPFJet_VBF_(0),
    HLT_L1_TripleJet_VBF_(0),
    HLT_QuadJet45_TripleBTagCSV_p087_(0),
    HLT_QuadJet45_DoubleBTagCSV_p087_(0),
    HLT_DoubleJet90_Double30_TripleBTagCSV_p087_(0),
    HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_(0),
    HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_(0),
    HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_(0),
    HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172_(0),
    HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_(0),
    HLT_DoubleJetsC100_SingleBTagCSV_p026_(0),
    HLT_DoubleJetsC100_SingleBTagCSV_p014_(0),
    HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350_(0),
    HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350_(0),
    HLT_Photon135_PFMET100_(0),
    HLT_Photon20_CaloIdVL_IsoL_(0),
    HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40_(0),
    HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_(0),
    HLT_Photon250_NoHE_(0),
    HLT_Photon300_NoHE_(0),
    HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60_(0),
    HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15_(0),
    HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_(0),
    HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_(0),
    HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_(0),
    HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_(0),
    HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40_(0),
    HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_(0),
    HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40_(0),
    HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_(0),
    HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_(0),
    HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_(0),
    HLT_Mu8_TrkIsoVVL_(0),
    HLT_Mu17_TrkIsoVVL_(0),
    HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    HLT_BTagMu_DiJet20_Mu5_(0),
    HLT_BTagMu_DiJet40_Mu5_(0),
    HLT_BTagMu_DiJet70_Mu5_(0),
    HLT_BTagMu_DiJet110_Mu5_(0),
    HLT_BTagMu_DiJet170_Mu5_(0),
    HLT_BTagMu_Jet300_Mu5_(0),
    HLT_BTagMu_AK8Jet300_Mu5_(0),
    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_(0),
    HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_(0),
    HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_(0),
    HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_(0),
    HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL_(0),
    HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL_(0),
    HLT_Mu8_DiEle12_CaloIdL_TrackIdL_(0),
    HLT_Mu12_Photon25_CaloIdL_(0),
    HLT_Mu12_Photon25_CaloIdL_L1ISO_(0),
    HLT_Mu12_Photon25_CaloIdL_L1OR_(0),
    HLT_Mu17_Photon22_CaloIdL_L1ISO_(0),
    HLT_Mu17_Photon30_CaloIdL_L1ISO_(0),
    HLT_Mu17_Photon35_CaloIdL_L1ISO_(0),
    HLT_DiMu9_Ele9_CaloIdL_TrackIdL_(0),
    HLT_TripleMu_5_3_3_(0),
    HLT_TripleMu_12_10_5_(0),
    HLT_Mu3er_PFHT140_PFMET125_(0),
    HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067_(0),
    HLT_Mu6_PFHT200_PFMET100_(0),
    HLT_Mu14er_PFMET100_(0),
    HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Ele17_CaloIdL_GsfTrkIdVL_(0),
    HLT_Ele17_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Ele23_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_(0),
    HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_(0),
    HLT_Photon22_(0),
    HLT_Photon30_(0),
    HLT_Photon36_(0),
    HLT_Photon50_(0),
    HLT_Photon75_(0),
    HLT_Photon90_(0),
    HLT_Photon120_(0),
    HLT_Photon175_(0),
    HLT_Photon165_HE10_(0),
    HLT_Photon22_R9Id90_HE10_IsoM_(0),
    HLT_Photon30_R9Id90_HE10_IsoM_(0),
    HLT_Photon36_R9Id90_HE10_IsoM_(0),
    HLT_Photon50_R9Id90_HE10_IsoM_(0),
    HLT_Photon75_R9Id90_HE10_IsoM_(0),
    HLT_Photon90_R9Id90_HE10_IsoM_(0),
    HLT_Photon120_R9Id90_HE10_IsoM_(0),
    HLT_Photon165_R9Id90_HE10_IsoM_(0),
    HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_(0),
    HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_(0),
    HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_(0),
    HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55_(0),
    HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_(0),
    HLT_Dimuon0_Jpsi_Muon_(0),
    HLT_Dimuon0_Upsilon_Muon_(0),
    HLT_QuadMuon0_Dimuon0_Jpsi_(0),
    HLT_QuadMuon0_Dimuon0_Upsilon_(0),
    HLT_Rsq0p25_Calo_(0),
    HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo_(0),
    HLT_RsqMR240_Rsq0p09_MR200_Calo_(0),
    HLT_Rsq0p25_(0),
    HLT_Rsq0p30_(0),
    HLT_RsqMR240_Rsq0p09_MR200_(0),
    HLT_RsqMR240_Rsq0p09_MR200_4jet_(0),
    HLT_RsqMR270_Rsq0p09_MR200_(0),
    HLT_RsqMR270_Rsq0p09_MR200_4jet_(0),
    HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200_(0),
    HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_(0),
    HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_(0),
    HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_(0),
    HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_(0),
    HLT_HT200_DisplacedDijet40_DisplacedTrack_(0),
    HLT_HT250_DisplacedDijet40_DisplacedTrack_(0),
    HLT_HT350_DisplacedDijet40_DisplacedTrack_(0),
    HLT_HT350_DisplacedDijet80_DisplacedTrack_(0),
    HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack_(0),
    HLT_HT350_DisplacedDijet40_Inclusive_(0),
    HLT_HT400_DisplacedDijet40_Inclusive_(0),
    HLT_HT500_DisplacedDijet40_Inclusive_(0),
    HLT_HT550_DisplacedDijet40_Inclusive_(0),
    HLT_HT550_DisplacedDijet80_Inclusive_(0),
    HLT_HT650_DisplacedDijet80_Inclusive_(0),
    HLT_HT750_DisplacedDijet80_Inclusive_(0),
    HLT_VBF_DisplacedJet40_DisplacedTrack_(0),
    HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5_(0),
    HLT_VBF_DisplacedJet40_TightID_DisplacedTrack_(0),
    HLT_VBF_DisplacedJet40_Hadronic_(0),
    HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack_(0),
    HLT_VBF_DisplacedJet40_TightID_Hadronic_(0),
    HLT_VBF_DisplacedJet40_VTightID_Hadronic_(0),
    HLT_VBF_DisplacedJet40_VVTightID_Hadronic_(0),
    HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_(0),
    HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_(0),
    HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_(0),
    HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_(0),
    HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_(0),
    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_(0),
    HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight_(0),
    HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight_(0),
    HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_(0),
    HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_(0),
    HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_(0),
    HLT_Photon90_CaloIdL_PFHT500_(0),
    HLT_DoubleMu8_Mass8_PFHT250_(0),
    HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250_(0),
    HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250_(0),
    HLT_DoubleMu8_Mass8_PFHT300_(0),
    HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_(0),
    HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_(0),
    HLT_Mu10_CentralPFJet30_BTagCSV_p13_(0),
    HLT_DoubleMu3_PFMET50_(0),
    HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_(0),
    HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_(0),
    HLT_Ele15_IsoVVVL_PFHT350_PFMET50_(0),
    HLT_Ele15_IsoVVVL_PFHT600_(0),
    HLT_Ele15_IsoVVVL_PFHT350_(0),
    HLT_Ele15_IsoVVVL_PFHT400_PFMET50_(0),
    HLT_Ele15_IsoVVVL_PFHT400_(0),
    HLT_Ele50_IsoVVVL_PFHT400_(0),
    HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_(0),
    HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_(0),
    HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_(0),
    HLT_Mu15_IsoVVVL_PFHT350_PFMET50_(0),
    HLT_Mu15_IsoVVVL_PFHT600_(0),
    HLT_Mu15_IsoVVVL_PFHT350_(0),
    HLT_Mu15_IsoVVVL_PFHT400_PFMET50_(0),
    HLT_Mu15_IsoVVVL_PFHT400_(0),
    HLT_Mu50_IsoVVVL_PFHT400_(0),
    HLT_Dimuon16_Jpsi_(0),
    HLT_Dimuon10_Jpsi_Barrel_(0),
    HLT_Dimuon8_PsiPrime_Barrel_(0),
    HLT_Dimuon8_Upsilon_Barrel_(0),
    HLT_Dimuon0_Phi_Barrel_(0),
    HLT_Mu16_TkMu0_dEta18_Onia_(0),
    HLT_Mu16_TkMu0_dEta18_Phi_(0),
    HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_(0),
    HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_(0),
    HLT_Mu8_(0),
    HLT_Mu17_(0),
    HLT_Mu3_PFJet40_(0),
    HLT_Ele8_CaloIdM_TrackIdM_PFJet30_(0),
    HLT_Ele12_CaloIdM_TrackIdM_PFJet30_(0),
    HLT_Ele17_CaloIdM_TrackIdM_PFJet30_(0),
    HLT_Ele23_CaloIdM_TrackIdM_PFJet30_(0),
    HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_(0),
    HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_(0),
    HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_(0),
    HLT_PFHT450_SixJet40_BTagCSV_p056_(0),
    HLT_PFHT400_SixJet30_(0),
    HLT_PFHT450_SixJet40_(0),
    HLT_Ele115_CaloIdVT_GsfTrkIdT_(0),
    HLT_Mu55_(0),
    HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_(0),
    HLT_Photon90_CaloIdL_PFHT600_(0),
    HLT_PixelTracks_Multiplicity60ForEndOfFill_(0),
    HLT_PixelTracks_Multiplicity85ForEndOfFill_(0),
    HLT_PixelTracks_Multiplicity110ForEndOfFill_(0),
    HLT_PixelTracks_Multiplicity135ForEndOfFill_(0),
    HLT_PixelTracks_Multiplicity160ForEndOfFill_(0),
    HLT_FullTracks_Multiplicity80_(0),
    HLT_FullTracks_Multiplicity100_(0),
    HLT_FullTracks_Multiplicity130_(0),
    HLT_FullTracks_Multiplicity150_(0),
    HLT_ECALHT800_(0),
    HLT_DiSC30_18_EIso_AND_HE_Mass70_(0),
    HLT_Photon125_(0),
    HLT_MET100_(0),
    HLT_MET150_(0),
    HLT_MET200_(0),
    HLT_Ele27_HighEta_Ele20_Mass55_(0),
    HLT_L1FatEvents_(0),
    HLT_Physics_(0),
    HLT_L1FatEvents_part0_(0),
    HLT_L1FatEvents_part1_(0),
    HLT_L1FatEvents_part2_(0),
    HLT_L1FatEvents_part3_(0),
    HLT_Random_(0),
    HLT_ZeroBias_(0),
    HLT_AK4CaloJet30_(0),
    HLT_AK4CaloJet40_(0),
    HLT_AK4CaloJet50_(0),
    HLT_AK4CaloJet80_(0),
    HLT_AK4CaloJet100_(0),
    HLT_AK4PFJet30_(0),
    HLT_AK4PFJet50_(0),
    HLT_AK4PFJet80_(0),
    HLT_AK4PFJet100_(0),
    HLT_HISinglePhoton10_(0),
    HLT_HISinglePhoton15_(0),
    HLT_HISinglePhoton20_(0),
    HLT_HISinglePhoton40_(0),
    HLT_HISinglePhoton60_(0),
    HLT_EcalCalibration_(0),
    HLT_HcalCalibration_(0),
    HLT_GlobalRunHPDNoise_(0),
    HLT_L1BptxMinus_(0),
    HLT_L1BptxPlus_(0),
    HLT_L1NotBptxOR_(0),
    HLT_L1BeamGasMinus_(0),
    HLT_L1BeamGasPlus_(0),
    HLT_L1BptxXOR_(0),
    HLT_L1MinimumBiasHF_OR_(0),
    HLT_L1MinimumBiasHF_AND_(0),
    HLT_HcalNZS_(0),
    HLT_HcalPhiSym_(0),
    HLT_HcalIsolatedbunch_(0),
    HLT_ZeroBias_FirstCollisionAfterAbortGap_(0),
    HLT_ZeroBias_FirstCollisionAfterAbortGap_copy_(0),
    HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS_(0),
    HLT_ZeroBias_IsolatedBunches_(0),
    HLT_ZeroBias_FirstCollisionInTrain_(0),
    HLT_ZeroBias_FirstBXAfterTrain_(0),
    HLT_Photon500_(0),
    HLT_Photon600_(0),
    HLT_Mu300_(0),
    HLT_Mu350_(0),
    HLT_MET250_(0),
    HLT_MET300_(0),
    HLT_MET600_(0),
    HLT_MET700_(0),
    HLT_PFMET300_(0),
    HLT_PFMET400_(0),
    HLT_PFMET500_(0),
    HLT_PFMET600_(0),
    HLT_Ele250_CaloIdVT_GsfTrkIdT_(0),
    HLT_Ele300_CaloIdVT_GsfTrkIdT_(0),
    HLT_HT2000_(0),
    HLT_HT2500_(0),
    HLT_IsoTrackHE_(0),
    HLT_IsoTrackHB_(0),
    HLTriggerFinalPath(0),
    Flag_HBHENoiseFilter_(0),
    Flag_HBHENoiseIsoFilter_(0),
    Flag_CSCTightHaloFilter_(0),
    Flag_CSCTightHaloTrkMuUnvetoFilter_(0),
    Flag_CSCTightHalo2015Filter_(0),
    Flag_globalTightHalo2016Filter_(0),
    Flag_globalSuperTightHalo2016Filter_(0),
    Flag_HcalStripHaloFilter_(0),
    Flag_hcalLaserEventFilter_(0),
    Flag_EcalDeadCellTriggerPrimitiveFilter_(0),
    Flag_EcalDeadCellBoundaryEnergyFilter_(0),
    Flag_goodVertices_(0),
    Flag_eeBadScFilter_(0),
    Flag_ecalLaserCorrFilter_(0),
    Flag_trkPOGFilters_(0),
    Flag_chargedHadronTrackResolutionFilter_(0),
    Flag_muonBadTrackFilter_(0),
    Flag_trkPOG_manystripclus53X_(0),
    Flag_trkPOG_toomanystripclus53X_(0),
    Flag_trkPOG_logErrorTooManyClusters_(0),
    Flag_METFilters_(0),
    are_Electron_loaded_(0), Electron_(),
    are_FatJet_loaded_(0), FatJet_(),
    are_GenJetAK8_loaded_(0), GenJetAK8_(),
    are_GenJet_loaded_(0), GenJet_(),
    are_GenPart_loaded_(0), GenPart_(),
    are_GenVisTau_loaded_(0), GenVisTau_(),
    are_LHEPdfWeight_loaded_(0), LHEPdfWeight_(),
    are_LHEScaleWeight_loaded_(0), LHEScaleWeight_(),
    are_Jet_loaded_(0), Jet_(),
    are_Muon_loaded_(0), Muon_(),
    are_Photon_loaded_(0), Photon_(),
    are_GenDressedLepton_loaded_(0), GenDressedLepton_(),
    are_SoftActivityJet_loaded_(0), SoftActivityJet_(),
    are_SubJet_loaded_(0), SubJet_(),
    are_Tau_loaded_(0), Tau_(),
    are_TrigObj_loaded_(0), TrigObj_(),
    are_OtherPV_loaded_(0), OtherPV_(),
    are_SV_loaded_(0), SV_()
  {
    tree_ = tree;
    current_entry_ = -1;
    entries_ = tree_->GetEntries();
    tree_->SetBranchStatus("*",0); 
    tree_->SetBranchStatus("run", 1); tree_->SetBranchAddress("run", &run);
    tree_->SetBranchStatus("luminosityBlock", 1); tree_->SetBranchAddress("luminosityBlock", &luminosityBlock);
    tree_->SetBranchStatus("event", 1); tree_->SetBranchAddress("event", &event);
    tree_->SetBranchStatus("nElectron", 1); tree_->SetBranchAddress("nElectron", &nElectron);
    tree_->SetBranchStatus("nFatJet", 1); tree_->SetBranchAddress("nFatJet", &nFatJet);
    tree_->SetBranchStatus("nGenJetAK8", 1); tree_->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
    tree_->SetBranchStatus("nGenJet", 1); tree_->SetBranchAddress("nGenJet", &nGenJet);
    tree_->SetBranchStatus("nGenPart", 1); tree_->SetBranchAddress("nGenPart", &nGenPart);
    tree_->SetBranchStatus("nGenVisTau", 1); tree_->SetBranchAddress("nGenVisTau", &nGenVisTau);
    tree_->SetBranchStatus("genWeight", 1); tree_->SetBranchAddress("genWeight", &genWeight);
    tree_->SetBranchStatus("nLHEPdfWeight", 1); tree_->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight);
    tree_->SetBranchStatus("LHEPdfWeight", 1); tree_->SetBranchAddress("LHEPdfWeight", &LHEPdfWeight);
    tree_->SetBranchStatus("nLHEScaleWeight", 1); tree_->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight);
    tree_->SetBranchStatus("LHEScaleWeight", 1); tree_->SetBranchAddress("LHEScaleWeight", &LHEScaleWeight);
    tree_->SetBranchStatus("nJet", 1); tree_->SetBranchAddress("nJet", &nJet);
    tree_->SetBranchStatus("nMuon", 1); tree_->SetBranchAddress("nMuon", &nMuon);
    tree_->SetBranchStatus("nPhoton", 1); tree_->SetBranchAddress("nPhoton", &nPhoton);
    tree_->SetBranchStatus("fixedGridRhoFastjetAll", 1); tree_->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll);
    tree_->SetBranchStatus("fixedGridRhoFastjetCentralCalo", 1); tree_->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo);
    tree_->SetBranchStatus("fixedGridRhoFastjetCentralNeutral", 1); tree_->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral);
    tree_->SetBranchStatus("nGenDressedLepton", 1); tree_->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton);
    tree_->SetBranchStatus("nSoftActivityJet", 1); tree_->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet);
    tree_->SetBranchStatus("SoftActivityJetHT", 1); tree_->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT);
    tree_->SetBranchStatus("SoftActivityJetHT10", 1); tree_->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10);
    tree_->SetBranchStatus("SoftActivityJetHT2", 1); tree_->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2);
    tree_->SetBranchStatus("SoftActivityJetHT5", 1); tree_->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5);
    tree_->SetBranchStatus("SoftActivityJetNjets10", 1); tree_->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10);
    tree_->SetBranchStatus("SoftActivityJetNjets2", 1); tree_->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2);
    tree_->SetBranchStatus("SoftActivityJetNjets5", 1); tree_->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5);
    tree_->SetBranchStatus("nSubJet", 1); tree_->SetBranchAddress("nSubJet", &nSubJet);
    tree_->SetBranchStatus("nTau", 1); tree_->SetBranchAddress("nTau", &nTau);
    tree_->SetBranchStatus("nTrigObj", 1); tree_->SetBranchAddress("nTrigObj", &nTrigObj);
    tree_->SetBranchStatus("genTtbarId", 1); tree_->SetBranchAddress("genTtbarId", &genTtbarId);
    tree_->SetBranchStatus("nOtherPV", 1); tree_->SetBranchAddress("nOtherPV", &nOtherPV);
    tree_->SetBranchStatus("nSV", 1); tree_->SetBranchAddress("nSV", &nSV);
    tree_->SetBranchStatus("HLTriggerFirstPath", 1); tree_->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath);
    tree_->SetBranchStatus("HLTriggerFinalPath", 1); tree_->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath);
  }

  ~URStreamer()
  {
    //{ EVT_DESTROY }
  }

  bool next(){
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    current_entry_++;
    if(current_entry_ < entries_){
      tree_->GetEntry(current_entry_);
      return true;
    }
    return false;
  }

  void loadElectron(){
    if(!are_Electron_loaded_){
      tree_->SetBranchStatus("Electron_deltaEtaSC", 1); tree_->SetBranchAddress("Electron_deltaEtaSC", &Electron_deltaEtaSC_);
      tree_->SetBranchStatus("Electron_dr03EcalRecHitSumEt", 1); tree_->SetBranchAddress("Electron_dr03EcalRecHitSumEt", &Electron_dr03EcalRecHitSumEt_);
      tree_->SetBranchStatus("Electron_dr03HcalDepth1TowerSumEt", 1); tree_->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", &Electron_dr03HcalDepth1TowerSumEt_);
      tree_->SetBranchStatus("Electron_dr03TkSumPt", 1); tree_->SetBranchAddress("Electron_dr03TkSumPt", &Electron_dr03TkSumPt_);
      tree_->SetBranchStatus("Electron_dxy", 1); tree_->SetBranchAddress("Electron_dxy", &Electron_dxy_);
      tree_->SetBranchStatus("Electron_dxyErr", 1); tree_->SetBranchAddress("Electron_dxyErr", &Electron_dxyErr_);
      tree_->SetBranchStatus("Electron_dz", 1); tree_->SetBranchAddress("Electron_dz", &Electron_dz_);
      tree_->SetBranchStatus("Electron_dzErr", 1); tree_->SetBranchAddress("Electron_dzErr", &Electron_dzErr_);
      tree_->SetBranchStatus("Electron_eCorr", 1); tree_->SetBranchAddress("Electron_eCorr", &Electron_eCorr_);
      tree_->SetBranchStatus("Electron_eInvMinusPInv", 1); tree_->SetBranchAddress("Electron_eInvMinusPInv", &Electron_eInvMinusPInv_);
      tree_->SetBranchStatus("Electron_energyErr", 1); tree_->SetBranchAddress("Electron_energyErr", &Electron_energyErr_);
      tree_->SetBranchStatus("Electron_eta", 1); tree_->SetBranchAddress("Electron_eta", &Electron_eta_);
      tree_->SetBranchStatus("Electron_hoe", 1); tree_->SetBranchAddress("Electron_hoe", &Electron_hoe_);
      tree_->SetBranchStatus("Electron_ip3d", 1); tree_->SetBranchAddress("Electron_ip3d", &Electron_ip3d_);
      tree_->SetBranchStatus("Electron_mass", 1); tree_->SetBranchAddress("Electron_mass", &Electron_mass_);
      tree_->SetBranchStatus("Electron_miniPFRelIso_all", 1); tree_->SetBranchAddress("Electron_miniPFRelIso_all", &Electron_miniPFRelIso_all_);
      tree_->SetBranchStatus("Electron_miniPFRelIso_chg", 1); tree_->SetBranchAddress("Electron_miniPFRelIso_chg", &Electron_miniPFRelIso_chg_);
      tree_->SetBranchStatus("Electron_mvaSpring16GP", 1); tree_->SetBranchAddress("Electron_mvaSpring16GP", &Electron_mvaSpring16GP_);
      tree_->SetBranchStatus("Electron_mvaSpring16HZZ", 1); tree_->SetBranchAddress("Electron_mvaSpring16HZZ", &Electron_mvaSpring16HZZ_);
      tree_->SetBranchStatus("Electron_pfRelIso03_all", 1); tree_->SetBranchAddress("Electron_pfRelIso03_all", &Electron_pfRelIso03_all_);
      tree_->SetBranchStatus("Electron_pfRelIso03_chg", 1); tree_->SetBranchAddress("Electron_pfRelIso03_chg", &Electron_pfRelIso03_chg_);
      tree_->SetBranchStatus("Electron_phi", 1); tree_->SetBranchAddress("Electron_phi", &Electron_phi_);
      tree_->SetBranchStatus("Electron_pt", 1); tree_->SetBranchAddress("Electron_pt", &Electron_pt_);
      tree_->SetBranchStatus("Electron_r9", 1); tree_->SetBranchAddress("Electron_r9", &Electron_r9_);
      tree_->SetBranchStatus("Electron_sieie", 1); tree_->SetBranchAddress("Electron_sieie", &Electron_sieie_);
      tree_->SetBranchStatus("Electron_sip3d", 1); tree_->SetBranchAddress("Electron_sip3d", &Electron_sip3d_);
      tree_->SetBranchStatus("Electron_mvaTTH", 1); tree_->SetBranchAddress("Electron_mvaTTH", &Electron_mvaTTH_);
      tree_->SetBranchStatus("Electron_charge", 1); tree_->SetBranchAddress("Electron_charge", &Electron_charge_);
      tree_->SetBranchStatus("Electron_cutBased", 1); tree_->SetBranchAddress("Electron_cutBased", &Electron_cutBased_);
      tree_->SetBranchStatus("Electron_cutBased_HLTPreSel", 1); tree_->SetBranchAddress("Electron_cutBased_HLTPreSel", &Electron_cutBased_HLTPreSel_);
      tree_->SetBranchStatus("Electron_jetIdx", 1); tree_->SetBranchAddress("Electron_jetIdx", &Electron_jetIdx_);
      tree_->SetBranchStatus("Electron_pdgId", 1); tree_->SetBranchAddress("Electron_pdgId", &Electron_pdgId_);
      tree_->SetBranchStatus("Electron_photonIdx", 1); tree_->SetBranchAddress("Electron_photonIdx", &Electron_photonIdx_);
      tree_->SetBranchStatus("Electron_tightCharge", 1); tree_->SetBranchAddress("Electron_tightCharge", &Electron_tightCharge_);
      tree_->SetBranchStatus("Electron_vidNestedWPBitmap", 1); tree_->SetBranchAddress("Electron_vidNestedWPBitmap", &Electron_vidNestedWPBitmap_);
      tree_->SetBranchStatus("Electron_convVeto", 1); tree_->SetBranchAddress("Electron_convVeto", &Electron_convVeto_);
      tree_->SetBranchStatus("Electron_cutBased_HEEP", 1); tree_->SetBranchAddress("Electron_cutBased_HEEP", &Electron_cutBased_HEEP_);
      tree_->SetBranchStatus("Electron_isPFcand", 1); tree_->SetBranchAddress("Electron_isPFcand", &Electron_isPFcand_);
      tree_->SetBranchStatus("Electron_lostHits", 1); tree_->SetBranchAddress("Electron_lostHits", &Electron_lostHits_);
      tree_->SetBranchStatus("Electron_mvaSpring16GP_WP80", 1); tree_->SetBranchAddress("Electron_mvaSpring16GP_WP80", &Electron_mvaSpring16GP_WP80_);
      tree_->SetBranchStatus("Electron_mvaSpring16GP_WP90", 1); tree_->SetBranchAddress("Electron_mvaSpring16GP_WP90", &Electron_mvaSpring16GP_WP90_);
      tree_->SetBranchStatus("Electron_mvaSpring16HZZ_WPL", 1); tree_->SetBranchAddress("Electron_mvaSpring16HZZ_WPL", &Electron_mvaSpring16HZZ_WPL_);
      tree_->SetBranchStatus("Electron_genPartIdx", 1); tree_->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx_);
      tree_->SetBranchStatus("Electron_genPartFlav", 1); tree_->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav_);
      tree_->SetBranchStatus("Electron_cleanmask", 1); tree_->SetBranchAddress("Electron_cleanmask", &Electron_cleanmask_);
      are_Electron_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadFatjet(){
    if(!are_FatJet_loaded_){
      tree_->SetBranchStatus("FatJet_area", 1); tree_->SetBranchAddress("FatJet_area", &FatJet_area_);
      tree_->SetBranchStatus("FatJet_btagCMVA", 1); tree_->SetBranchAddress("FatJet_btagCMVA", &FatJet_btagCMVA_);
      tree_->SetBranchStatus("FatJet_btagCSVV2", 1); tree_->SetBranchAddress("FatJet_btagCSVV2", &FatJet_btagCSVV2_);
      tree_->SetBranchStatus("FatJet_btagDeepB", 1); tree_->SetBranchAddress("FatJet_btagDeepB", &FatJet_btagDeepB_);
      tree_->SetBranchStatus("FatJet_btagHbb", 1); tree_->SetBranchAddress("FatJet_btagHbb", &FatJet_btagHbb_);
      tree_->SetBranchStatus("FatJet_eta", 1); tree_->SetBranchAddress("FatJet_eta", &FatJet_eta_);
      tree_->SetBranchStatus("FatJet_mass", 1); tree_->SetBranchAddress("FatJet_mass", &FatJet_mass_);
      tree_->SetBranchStatus("FatJet_msoftdrop", 1); tree_->SetBranchAddress("FatJet_msoftdrop", &FatJet_msoftdrop_);
      tree_->SetBranchStatus("FatJet_msoftdrop_chs", 1); tree_->SetBranchAddress("FatJet_msoftdrop_chs", &FatJet_msoftdrop_chs_);
      tree_->SetBranchStatus("FatJet_n2b1", 1); tree_->SetBranchAddress("FatJet_n2b1", &FatJet_n2b1_);
      tree_->SetBranchStatus("FatJet_n3b1", 1); tree_->SetBranchAddress("FatJet_n3b1", &FatJet_n3b1_);
      tree_->SetBranchStatus("FatJet_phi", 1); tree_->SetBranchAddress("FatJet_phi", &FatJet_phi_);
      tree_->SetBranchStatus("FatJet_pt", 1); tree_->SetBranchAddress("FatJet_pt", &FatJet_pt_);
      tree_->SetBranchStatus("FatJet_tau1", 1); tree_->SetBranchAddress("FatJet_tau1", &FatJet_tau1_);
      tree_->SetBranchStatus("FatJet_tau2", 1); tree_->SetBranchAddress("FatJet_tau2", &FatJet_tau2_);
      tree_->SetBranchStatus("FatJet_tau3", 1); tree_->SetBranchAddress("FatJet_tau3", &FatJet_tau3_);
      tree_->SetBranchStatus("FatJet_tau4", 1); tree_->SetBranchAddress("FatJet_tau4", &FatJet_tau4_);
      tree_->SetBranchStatus("FatJet_jetId", 1); tree_->SetBranchAddress("FatJet_jetId", &FatJet_jetId_);
      tree_->SetBranchStatus("FatJet_subJetIdx1", 1); tree_->SetBranchAddress("FatJet_subJetIdx1", &FatJet_subJetIdx1_);
      tree_->SetBranchStatus("FatJet_subJetIdx2", 1); tree_->SetBranchAddress("FatJet_subJetIdx2", &FatJet_subJetIdx2_);
      are_FatJet_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenjetak8(){
    if(!are_GenJetAK8_loaded_){
      tree_->SetBranchStatus("GenJetAK8_eta", 1); tree_->SetBranchAddress("GenJetAK8_eta", &GenJetAK8_eta_);
      tree_->SetBranchStatus("GenJetAK8_mass", 1); tree_->SetBranchAddress("GenJetAK8_mass", &GenJetAK8_mass_);
      tree_->SetBranchStatus("GenJetAK8_phi", 1); tree_->SetBranchAddress("GenJetAK8_phi", &GenJetAK8_phi_);
      tree_->SetBranchStatus("GenJetAK8_pt", 1); tree_->SetBranchAddress("GenJetAK8_pt", &GenJetAK8_pt_);
      tree_->SetBranchStatus("GenJetAK8_partonFlavour", 1); tree_->SetBranchAddress("GenJetAK8_partonFlavour", &GenJetAK8_partonFlavour_);
      tree_->SetBranchStatus("GenJetAK8_hadronFlavour", 1); tree_->SetBranchAddress("GenJetAK8_hadronFlavour", &GenJetAK8_hadronFlavour_);
      are_GenJetAK8_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenjet(){
    if(!are_GenJet_loaded_){
      tree_->SetBranchStatus("GenJetAK8_eta", 1); tree_->SetBranchAddress("GenJetAK8_eta", &GenJetAK8_eta_);
      tree_->SetBranchStatus("GenJetAK8_mass", 1); tree_->SetBranchAddress("GenJetAK8_mass", &GenJetAK8_mass_);
      tree_->SetBranchStatus("GenJetAK8_phi", 1); tree_->SetBranchAddress("GenJetAK8_phi", &GenJetAK8_phi_);
      tree_->SetBranchStatus("GenJetAK8_pt", 1); tree_->SetBranchAddress("GenJetAK8_pt", &GenJetAK8_pt_);
      tree_->SetBranchStatus("GenJet_eta", 1); tree_->SetBranchAddress("GenJet_eta", &GenJet_eta_);
      tree_->SetBranchStatus("GenJet_mass", 1); tree_->SetBranchAddress("GenJet_mass", &GenJet_mass_);
      tree_->SetBranchStatus("GenJet_phi", 1); tree_->SetBranchAddress("GenJet_phi", &GenJet_phi_);
      tree_->SetBranchStatus("GenJet_pt", 1); tree_->SetBranchAddress("GenJet_pt", &GenJet_pt_);
      tree_->SetBranchStatus("GenJetAK8_partonFlavour", 1); tree_->SetBranchAddress("GenJetAK8_partonFlavour", &GenJetAK8_partonFlavour_);
      tree_->SetBranchStatus("GenJetAK8_hadronFlavour", 1); tree_->SetBranchAddress("GenJetAK8_hadronFlavour", &GenJetAK8_hadronFlavour_);
      tree_->SetBranchStatus("GenJet_partonFlavour", 1); tree_->SetBranchAddress("GenJet_partonFlavour", &GenJet_partonFlavour_);
      tree_->SetBranchStatus("GenJet_hadronFlavour", 1); tree_->SetBranchAddress("GenJet_hadronFlavour", &GenJet_hadronFlavour_);
      are_GenJet_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenpart(){
    if(!are_GenPart_loaded_){
      tree_->SetBranchStatus("GenPart_eta", 1); tree_->SetBranchAddress("GenPart_eta", &GenPart_eta_);
      tree_->SetBranchStatus("GenPart_mass", 1); tree_->SetBranchAddress("GenPart_mass", &GenPart_mass_);
      tree_->SetBranchStatus("GenPart_phi", 1); tree_->SetBranchAddress("GenPart_phi", &GenPart_phi_);
      tree_->SetBranchStatus("GenPart_pt", 1); tree_->SetBranchAddress("GenPart_pt", &GenPart_pt_);
      tree_->SetBranchStatus("GenPart_genPartIdxMother", 1); tree_->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother_);
      tree_->SetBranchStatus("GenPart_pdgId", 1); tree_->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId_);
      tree_->SetBranchStatus("GenPart_status", 1); tree_->SetBranchAddress("GenPart_status", &GenPart_status_);
      tree_->SetBranchStatus("GenPart_statusFlags", 1); tree_->SetBranchAddress("GenPart_statusFlags", &GenPart_statusFlags_);
      are_GenPart_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenvistau(){
    if(!are_GenVisTau_loaded_){
      tree_->SetBranchStatus("GenVisTau_eta", 1); tree_->SetBranchAddress("GenVisTau_eta", &GenVisTau_eta_);
      tree_->SetBranchStatus("GenVisTau_mass", 1); tree_->SetBranchAddress("GenVisTau_mass", &GenVisTau_mass_);
      tree_->SetBranchStatus("GenVisTau_phi", 1); tree_->SetBranchAddress("GenVisTau_phi", &GenVisTau_phi_);
      tree_->SetBranchStatus("GenVisTau_pt", 1); tree_->SetBranchAddress("GenVisTau_pt", &GenVisTau_pt_);
      tree_->SetBranchStatus("GenVisTau_charge", 1); tree_->SetBranchAddress("GenVisTau_charge", &GenVisTau_charge_);
      tree_->SetBranchStatus("GenVisTau_genPartIdxMother", 1); tree_->SetBranchAddress("GenVisTau_genPartIdxMother", &GenVisTau_genPartIdxMother_);
      tree_->SetBranchStatus("GenVisTau_status", 1); tree_->SetBranchAddress("GenVisTau_status", &GenVisTau_status_);
      are_GenVisTau_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadLhepdfweight(){
    if(!are_LHEPdfWeight_loaded_){
      tree_->SetBranchStatus("LHEPdfWeight", 1); tree_->SetBranchAddress("LHEPdfWeight", &LHEPdfWeight);
      are_LHEPdfWeight_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadLhescaleweight(){
    if(!are_LHEScaleWeight_loaded_){
      tree_->SetBranchStatus("LHEScaleWeight", 1); tree_->SetBranchAddress("LHEScaleWeight", &LHEScaleWeight);
      are_LHEScaleWeight_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadJet(){
    if(!are_Jet_loaded_){
      tree_->SetBranchStatus("Jet_area", 1); tree_->SetBranchAddress("Jet_area", &Jet_area_);
      tree_->SetBranchStatus("Jet_btagCMVA", 1); tree_->SetBranchAddress("Jet_btagCMVA", &Jet_btagCMVA_);
      tree_->SetBranchStatus("Jet_btagCSVV2", 1); tree_->SetBranchAddress("Jet_btagCSVV2", &Jet_btagCSVV2_);
      tree_->SetBranchStatus("Jet_btagDeepB", 1); tree_->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB_);
      tree_->SetBranchStatus("Jet_btagDeepC", 1); tree_->SetBranchAddress("Jet_btagDeepC", &Jet_btagDeepC_);
      tree_->SetBranchStatus("Jet_chEmEF", 1); tree_->SetBranchAddress("Jet_chEmEF", &Jet_chEmEF_);
      tree_->SetBranchStatus("Jet_chHEF", 1); tree_->SetBranchAddress("Jet_chHEF", &Jet_chHEF_);
      tree_->SetBranchStatus("Jet_eta", 1); tree_->SetBranchAddress("Jet_eta", &Jet_eta_);
      tree_->SetBranchStatus("Jet_mass", 1); tree_->SetBranchAddress("Jet_mass", &Jet_mass_);
      tree_->SetBranchStatus("Jet_neEmEF", 1); tree_->SetBranchAddress("Jet_neEmEF", &Jet_neEmEF_);
      tree_->SetBranchStatus("Jet_neHEF", 1); tree_->SetBranchAddress("Jet_neHEF", &Jet_neHEF_);
      tree_->SetBranchStatus("Jet_phi", 1); tree_->SetBranchAddress("Jet_phi", &Jet_phi_);
      tree_->SetBranchStatus("Jet_pt", 1); tree_->SetBranchAddress("Jet_pt", &Jet_pt_);
      tree_->SetBranchStatus("Jet_qgl", 1); tree_->SetBranchAddress("Jet_qgl", &Jet_qgl_);
      tree_->SetBranchStatus("Jet_rawFactor", 1); tree_->SetBranchAddress("Jet_rawFactor", &Jet_rawFactor_);
      tree_->SetBranchStatus("Jet_bReg", 1); tree_->SetBranchAddress("Jet_bReg", &Jet_bReg_);
      tree_->SetBranchStatus("Jet_electronIdx1", 1); tree_->SetBranchAddress("Jet_electronIdx1", &Jet_electronIdx1_);
      tree_->SetBranchStatus("Jet_electronIdx2", 1); tree_->SetBranchAddress("Jet_electronIdx2", &Jet_electronIdx2_);
      tree_->SetBranchStatus("Jet_jetId", 1); tree_->SetBranchAddress("Jet_jetId", &Jet_jetId_);
      tree_->SetBranchStatus("Jet_muonIdx1", 1); tree_->SetBranchAddress("Jet_muonIdx1", &Jet_muonIdx1_);
      tree_->SetBranchStatus("Jet_muonIdx2", 1); tree_->SetBranchAddress("Jet_muonIdx2", &Jet_muonIdx2_);
      tree_->SetBranchStatus("Jet_nConstituents", 1); tree_->SetBranchAddress("Jet_nConstituents", &Jet_nConstituents_);
      tree_->SetBranchStatus("Jet_nElectrons", 1); tree_->SetBranchAddress("Jet_nElectrons", &Jet_nElectrons_);
      tree_->SetBranchStatus("Jet_nMuons", 1); tree_->SetBranchAddress("Jet_nMuons", &Jet_nMuons_);
      tree_->SetBranchStatus("Jet_puId", 1); tree_->SetBranchAddress("Jet_puId", &Jet_puId_);
      tree_->SetBranchStatus("Jet_genJetIdx", 1); tree_->SetBranchAddress("Jet_genJetIdx", &Jet_genJetIdx_);
      tree_->SetBranchStatus("Jet_hadronFlavour", 1); tree_->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour_);
      tree_->SetBranchStatus("Jet_partonFlavour", 1); tree_->SetBranchAddress("Jet_partonFlavour", &Jet_partonFlavour_);
      tree_->SetBranchStatus("Jet_cleanmask", 1); tree_->SetBranchAddress("Jet_cleanmask", &Jet_cleanmask_);
      are_Jet_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadMuon(){
    if(!are_Muon_loaded_){
      tree_->SetBranchStatus("Muon_dxy", 1); tree_->SetBranchAddress("Muon_dxy", &Muon_dxy_);
      tree_->SetBranchStatus("Muon_dxyErr", 1); tree_->SetBranchAddress("Muon_dxyErr", &Muon_dxyErr_);
      tree_->SetBranchStatus("Muon_dz", 1); tree_->SetBranchAddress("Muon_dz", &Muon_dz_);
      tree_->SetBranchStatus("Muon_dzErr", 1); tree_->SetBranchAddress("Muon_dzErr", &Muon_dzErr_);
      tree_->SetBranchStatus("Muon_eta", 1); tree_->SetBranchAddress("Muon_eta", &Muon_eta_);
      tree_->SetBranchStatus("Muon_ip3d", 1); tree_->SetBranchAddress("Muon_ip3d", &Muon_ip3d_);
      tree_->SetBranchStatus("Muon_mass", 1); tree_->SetBranchAddress("Muon_mass", &Muon_mass_);
      tree_->SetBranchStatus("Muon_miniPFRelIso_all", 1); tree_->SetBranchAddress("Muon_miniPFRelIso_all", &Muon_miniPFRelIso_all_);
      tree_->SetBranchStatus("Muon_miniPFRelIso_chg", 1); tree_->SetBranchAddress("Muon_miniPFRelIso_chg", &Muon_miniPFRelIso_chg_);
      tree_->SetBranchStatus("Muon_pfRelIso03_all", 1); tree_->SetBranchAddress("Muon_pfRelIso03_all", &Muon_pfRelIso03_all_);
      tree_->SetBranchStatus("Muon_pfRelIso03_chg", 1); tree_->SetBranchAddress("Muon_pfRelIso03_chg", &Muon_pfRelIso03_chg_);
      tree_->SetBranchStatus("Muon_pfRelIso04_all", 1); tree_->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all_);
      tree_->SetBranchStatus("Muon_phi", 1); tree_->SetBranchAddress("Muon_phi", &Muon_phi_);
      tree_->SetBranchStatus("Muon_pt", 1); tree_->SetBranchAddress("Muon_pt", &Muon_pt_);
      tree_->SetBranchStatus("Muon_ptErr", 1); tree_->SetBranchAddress("Muon_ptErr", &Muon_ptErr_);
      tree_->SetBranchStatus("Muon_segmentComp", 1); tree_->SetBranchAddress("Muon_segmentComp", &Muon_segmentComp_);
      tree_->SetBranchStatus("Muon_sip3d", 1); tree_->SetBranchAddress("Muon_sip3d", &Muon_sip3d_);
      tree_->SetBranchStatus("Muon_mvaTTH", 1); tree_->SetBranchAddress("Muon_mvaTTH", &Muon_mvaTTH_);
      tree_->SetBranchStatus("Muon_charge", 1); tree_->SetBranchAddress("Muon_charge", &Muon_charge_);
      tree_->SetBranchStatus("Muon_jetIdx", 1); tree_->SetBranchAddress("Muon_jetIdx", &Muon_jetIdx_);
      tree_->SetBranchStatus("Muon_nStations", 1); tree_->SetBranchAddress("Muon_nStations", &Muon_nStations_);
      tree_->SetBranchStatus("Muon_nTrackerLayers", 1); tree_->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers_);
      tree_->SetBranchStatus("Muon_pdgId", 1); tree_->SetBranchAddress("Muon_pdgId", &Muon_pdgId_);
      tree_->SetBranchStatus("Muon_tightCharge", 1); tree_->SetBranchAddress("Muon_tightCharge", &Muon_tightCharge_);
      tree_->SetBranchStatus("Muon_highPtId", 1); tree_->SetBranchAddress("Muon_highPtId", &Muon_highPtId_);
      tree_->SetBranchStatus("Muon_isPFcand", 1); tree_->SetBranchAddress("Muon_isPFcand", &Muon_isPFcand_);
      tree_->SetBranchStatus("Muon_mediumId", 1); tree_->SetBranchAddress("Muon_mediumId", &Muon_mediumId_);
      tree_->SetBranchStatus("Muon_softId", 1); tree_->SetBranchAddress("Muon_softId", &Muon_softId_);
      tree_->SetBranchStatus("Muon_tightId", 1); tree_->SetBranchAddress("Muon_tightId", &Muon_tightId_);
      tree_->SetBranchStatus("Muon_genPartIdx", 1); tree_->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx_);
      tree_->SetBranchStatus("Muon_genPartFlav", 1); tree_->SetBranchAddress("Muon_genPartFlav", &Muon_genPartFlav_);
      tree_->SetBranchStatus("Muon_cleanmask", 1); tree_->SetBranchAddress("Muon_cleanmask", &Muon_cleanmask_);
      are_Muon_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadPhoton(){
    if(!are_Photon_loaded_){
      tree_->SetBranchStatus("Photon_eCorr", 1); tree_->SetBranchAddress("Photon_eCorr", &Photon_eCorr_);
      tree_->SetBranchStatus("Photon_energyErr", 1); tree_->SetBranchAddress("Photon_energyErr", &Photon_energyErr_);
      tree_->SetBranchStatus("Photon_eta", 1); tree_->SetBranchAddress("Photon_eta", &Photon_eta_);
      tree_->SetBranchStatus("Photon_hoe", 1); tree_->SetBranchAddress("Photon_hoe", &Photon_hoe_);
      tree_->SetBranchStatus("Photon_mass", 1); tree_->SetBranchAddress("Photon_mass", &Photon_mass_);
      tree_->SetBranchStatus("Photon_mvaID", 1); tree_->SetBranchAddress("Photon_mvaID", &Photon_mvaID_);
      tree_->SetBranchStatus("Photon_pfRelIso03_all", 1); tree_->SetBranchAddress("Photon_pfRelIso03_all", &Photon_pfRelIso03_all_);
      tree_->SetBranchStatus("Photon_pfRelIso03_chg", 1); tree_->SetBranchAddress("Photon_pfRelIso03_chg", &Photon_pfRelIso03_chg_);
      tree_->SetBranchStatus("Photon_phi", 1); tree_->SetBranchAddress("Photon_phi", &Photon_phi_);
      tree_->SetBranchStatus("Photon_pt", 1); tree_->SetBranchAddress("Photon_pt", &Photon_pt_);
      tree_->SetBranchStatus("Photon_r9", 1); tree_->SetBranchAddress("Photon_r9", &Photon_r9_);
      tree_->SetBranchStatus("Photon_sieie", 1); tree_->SetBranchAddress("Photon_sieie", &Photon_sieie_);
      tree_->SetBranchStatus("Photon_charge", 1); tree_->SetBranchAddress("Photon_charge", &Photon_charge_);
      tree_->SetBranchStatus("Photon_cutBased", 1); tree_->SetBranchAddress("Photon_cutBased", &Photon_cutBased_);
      tree_->SetBranchStatus("Photon_electronIdx", 1); tree_->SetBranchAddress("Photon_electronIdx", &Photon_electronIdx_);
      tree_->SetBranchStatus("Photon_jetIdx", 1); tree_->SetBranchAddress("Photon_jetIdx", &Photon_jetIdx_);
      tree_->SetBranchStatus("Photon_pdgId", 1); tree_->SetBranchAddress("Photon_pdgId", &Photon_pdgId_);
      tree_->SetBranchStatus("Photon_vidNestedWPBitmap", 1); tree_->SetBranchAddress("Photon_vidNestedWPBitmap", &Photon_vidNestedWPBitmap_);
      tree_->SetBranchStatus("Photon_electronVeto", 1); tree_->SetBranchAddress("Photon_electronVeto", &Photon_electronVeto_);
      tree_->SetBranchStatus("Photon_mvaID_WP80", 1); tree_->SetBranchAddress("Photon_mvaID_WP80", &Photon_mvaID_WP80_);
      tree_->SetBranchStatus("Photon_mvaID_WP90", 1); tree_->SetBranchAddress("Photon_mvaID_WP90", &Photon_mvaID_WP90_);
      tree_->SetBranchStatus("Photon_pixelSeed", 1); tree_->SetBranchAddress("Photon_pixelSeed", &Photon_pixelSeed_);
      tree_->SetBranchStatus("Photon_genPartIdx", 1); tree_->SetBranchAddress("Photon_genPartIdx", &Photon_genPartIdx_);
      tree_->SetBranchStatus("Photon_genPartFlav", 1); tree_->SetBranchAddress("Photon_genPartFlav", &Photon_genPartFlav_);
      tree_->SetBranchStatus("Photon_cleanmask", 1); tree_->SetBranchAddress("Photon_cleanmask", &Photon_cleanmask_);
      are_Photon_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGendressedlepton(){
    if(!are_GenDressedLepton_loaded_){
      tree_->SetBranchStatus("GenDressedLepton_eta", 1); tree_->SetBranchAddress("GenDressedLepton_eta", &GenDressedLepton_eta_);
      tree_->SetBranchStatus("GenDressedLepton_mass", 1); tree_->SetBranchAddress("GenDressedLepton_mass", &GenDressedLepton_mass_);
      tree_->SetBranchStatus("GenDressedLepton_phi", 1); tree_->SetBranchAddress("GenDressedLepton_phi", &GenDressedLepton_phi_);
      tree_->SetBranchStatus("GenDressedLepton_pt", 1); tree_->SetBranchAddress("GenDressedLepton_pt", &GenDressedLepton_pt_);
      tree_->SetBranchStatus("GenDressedLepton_pdgId", 1); tree_->SetBranchAddress("GenDressedLepton_pdgId", &GenDressedLepton_pdgId_);
      are_GenDressedLepton_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadSoftactivityjet(){
    if(!are_SoftActivityJet_loaded_){
      tree_->SetBranchStatus("SoftActivityJet_eta", 1); tree_->SetBranchAddress("SoftActivityJet_eta", &SoftActivityJet_eta_);
      tree_->SetBranchStatus("SoftActivityJet_phi", 1); tree_->SetBranchAddress("SoftActivityJet_phi", &SoftActivityJet_phi_);
      tree_->SetBranchStatus("SoftActivityJet_pt", 1); tree_->SetBranchAddress("SoftActivityJet_pt", &SoftActivityJet_pt_);
      tree_->SetBranchStatus("SoftActivityJetHT", 1); tree_->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT);
      tree_->SetBranchStatus("SoftActivityJetHT10", 1); tree_->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10);
      tree_->SetBranchStatus("SoftActivityJetHT2", 1); tree_->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2);
      tree_->SetBranchStatus("SoftActivityJetHT5", 1); tree_->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5);
      tree_->SetBranchStatus("SoftActivityJetNjets10", 1); tree_->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10);
      tree_->SetBranchStatus("SoftActivityJetNjets2", 1); tree_->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2);
      tree_->SetBranchStatus("SoftActivityJetNjets5", 1); tree_->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5);
      are_SoftActivityJet_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadSubjet(){
    if(!are_SubJet_loaded_){
      tree_->SetBranchStatus("SubJet_btagCMVA", 1); tree_->SetBranchAddress("SubJet_btagCMVA", &SubJet_btagCMVA_);
      tree_->SetBranchStatus("SubJet_btagCSVV2", 1); tree_->SetBranchAddress("SubJet_btagCSVV2", &SubJet_btagCSVV2_);
      tree_->SetBranchStatus("SubJet_btagDeepB", 1); tree_->SetBranchAddress("SubJet_btagDeepB", &SubJet_btagDeepB_);
      tree_->SetBranchStatus("SubJet_eta", 1); tree_->SetBranchAddress("SubJet_eta", &SubJet_eta_);
      tree_->SetBranchStatus("SubJet_mass", 1); tree_->SetBranchAddress("SubJet_mass", &SubJet_mass_);
      tree_->SetBranchStatus("SubJet_n2b1", 1); tree_->SetBranchAddress("SubJet_n2b1", &SubJet_n2b1_);
      tree_->SetBranchStatus("SubJet_n3b1", 1); tree_->SetBranchAddress("SubJet_n3b1", &SubJet_n3b1_);
      tree_->SetBranchStatus("SubJet_phi", 1); tree_->SetBranchAddress("SubJet_phi", &SubJet_phi_);
      tree_->SetBranchStatus("SubJet_pt", 1); tree_->SetBranchAddress("SubJet_pt", &SubJet_pt_);
      tree_->SetBranchStatus("SubJet_tau1", 1); tree_->SetBranchAddress("SubJet_tau1", &SubJet_tau1_);
      tree_->SetBranchStatus("SubJet_tau2", 1); tree_->SetBranchAddress("SubJet_tau2", &SubJet_tau2_);
      tree_->SetBranchStatus("SubJet_tau3", 1); tree_->SetBranchAddress("SubJet_tau3", &SubJet_tau3_);
      tree_->SetBranchStatus("SubJet_tau4", 1); tree_->SetBranchAddress("SubJet_tau4", &SubJet_tau4_);
      are_SubJet_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadTau(){
    if(!are_Tau_loaded_){
      tree_->SetBranchStatus("Tau_chargedIso", 1); tree_->SetBranchAddress("Tau_chargedIso", &Tau_chargedIso_);
      tree_->SetBranchStatus("Tau_dxy", 1); tree_->SetBranchAddress("Tau_dxy", &Tau_dxy_);
      tree_->SetBranchStatus("Tau_dz", 1); tree_->SetBranchAddress("Tau_dz", &Tau_dz_);
      tree_->SetBranchStatus("Tau_eta", 1); tree_->SetBranchAddress("Tau_eta", &Tau_eta_);
      tree_->SetBranchStatus("Tau_footprintCorr", 1); tree_->SetBranchAddress("Tau_footprintCorr", &Tau_footprintCorr_);
      tree_->SetBranchStatus("Tau_leadTkDeltaEta", 1); tree_->SetBranchAddress("Tau_leadTkDeltaEta", &Tau_leadTkDeltaEta_);
      tree_->SetBranchStatus("Tau_leadTkDeltaPhi", 1); tree_->SetBranchAddress("Tau_leadTkDeltaPhi", &Tau_leadTkDeltaPhi_);
      tree_->SetBranchStatus("Tau_leadTkPtOverTauPt", 1); tree_->SetBranchAddress("Tau_leadTkPtOverTauPt", &Tau_leadTkPtOverTauPt_);
      tree_->SetBranchStatus("Tau_mass", 1); tree_->SetBranchAddress("Tau_mass", &Tau_mass_);
      tree_->SetBranchStatus("Tau_neutralIso", 1); tree_->SetBranchAddress("Tau_neutralIso", &Tau_neutralIso_);
      tree_->SetBranchStatus("Tau_phi", 1); tree_->SetBranchAddress("Tau_phi", &Tau_phi_);
      tree_->SetBranchStatus("Tau_photonsOutsideSignalCone", 1); tree_->SetBranchAddress("Tau_photonsOutsideSignalCone", &Tau_photonsOutsideSignalCone_);
      tree_->SetBranchStatus("Tau_pt", 1); tree_->SetBranchAddress("Tau_pt", &Tau_pt_);
      tree_->SetBranchStatus("Tau_puCorr", 1); tree_->SetBranchAddress("Tau_puCorr", &Tau_puCorr_);
      tree_->SetBranchStatus("Tau_rawAntiEle", 1); tree_->SetBranchAddress("Tau_rawAntiEle", &Tau_rawAntiEle_);
      tree_->SetBranchStatus("Tau_rawIso", 1); tree_->SetBranchAddress("Tau_rawIso", &Tau_rawIso_);
      tree_->SetBranchStatus("Tau_rawMVAnewDM", 1); tree_->SetBranchAddress("Tau_rawMVAnewDM", &Tau_rawMVAnewDM_);
      tree_->SetBranchStatus("Tau_rawMVAoldDM", 1); tree_->SetBranchAddress("Tau_rawMVAoldDM", &Tau_rawMVAoldDM_);
      tree_->SetBranchStatus("Tau_rawMVAoldDMdR03", 1); tree_->SetBranchAddress("Tau_rawMVAoldDMdR03", &Tau_rawMVAoldDMdR03_);
      tree_->SetBranchStatus("Tau_charge", 1); tree_->SetBranchAddress("Tau_charge", &Tau_charge_);
      tree_->SetBranchStatus("Tau_decayMode", 1); tree_->SetBranchAddress("Tau_decayMode", &Tau_decayMode_);
      tree_->SetBranchStatus("Tau_jetIdx", 1); tree_->SetBranchAddress("Tau_jetIdx", &Tau_jetIdx_);
      tree_->SetBranchStatus("Tau_rawAntiEleCat", 1); tree_->SetBranchAddress("Tau_rawAntiEleCat", &Tau_rawAntiEleCat_);
      tree_->SetBranchStatus("Tau_idAntiEle", 1); tree_->SetBranchAddress("Tau_idAntiEle", &Tau_idAntiEle_);
      tree_->SetBranchStatus("Tau_idAntiMu", 1); tree_->SetBranchAddress("Tau_idAntiMu", &Tau_idAntiMu_);
      tree_->SetBranchStatus("Tau_idDecayMode", 1); tree_->SetBranchAddress("Tau_idDecayMode", &Tau_idDecayMode_);
      tree_->SetBranchStatus("Tau_idDecayModeNewDMs", 1); tree_->SetBranchAddress("Tau_idDecayModeNewDMs", &Tau_idDecayModeNewDMs_);
      tree_->SetBranchStatus("Tau_idMVAnewDM", 1); tree_->SetBranchAddress("Tau_idMVAnewDM", &Tau_idMVAnewDM_);
      tree_->SetBranchStatus("Tau_idMVAoldDM", 1); tree_->SetBranchAddress("Tau_idMVAoldDM", &Tau_idMVAoldDM_);
      tree_->SetBranchStatus("Tau_idMVAoldDMdR03", 1); tree_->SetBranchAddress("Tau_idMVAoldDMdR03", &Tau_idMVAoldDMdR03_);
      tree_->SetBranchStatus("Tau_cleanmask", 1); tree_->SetBranchAddress("Tau_cleanmask", &Tau_cleanmask_);
      tree_->SetBranchStatus("Tau_genPartIdx", 1); tree_->SetBranchAddress("Tau_genPartIdx", &Tau_genPartIdx_);
      tree_->SetBranchStatus("Tau_genPartFlav", 1); tree_->SetBranchAddress("Tau_genPartFlav", &Tau_genPartFlav_);
      are_Tau_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadTrigobj(){
    if(!are_TrigObj_loaded_){
      tree_->SetBranchStatus("TrigObj_pt", 1); tree_->SetBranchAddress("TrigObj_pt", &TrigObj_pt_);
      tree_->SetBranchStatus("TrigObj_eta", 1); tree_->SetBranchAddress("TrigObj_eta", &TrigObj_eta_);
      tree_->SetBranchStatus("TrigObj_phi", 1); tree_->SetBranchAddress("TrigObj_phi", &TrigObj_phi_);
      tree_->SetBranchStatus("TrigObj_l1pt", 1); tree_->SetBranchAddress("TrigObj_l1pt", &TrigObj_l1pt_);
      tree_->SetBranchStatus("TrigObj_l1pt_2", 1); tree_->SetBranchAddress("TrigObj_l1pt_2", &TrigObj_l1pt_2_);
      tree_->SetBranchStatus("TrigObj_l2pt", 1); tree_->SetBranchAddress("TrigObj_l2pt", &TrigObj_l2pt_);
      tree_->SetBranchStatus("TrigObj_id", 1); tree_->SetBranchAddress("TrigObj_id", &TrigObj_id_);
      tree_->SetBranchStatus("TrigObj_l1iso", 1); tree_->SetBranchAddress("TrigObj_l1iso", &TrigObj_l1iso_);
      tree_->SetBranchStatus("TrigObj_l1charge", 1); tree_->SetBranchAddress("TrigObj_l1charge", &TrigObj_l1charge_);
      tree_->SetBranchStatus("TrigObj_filterBits", 1); tree_->SetBranchAddress("TrigObj_filterBits", &TrigObj_filterBits_);
      are_TrigObj_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadOtherpv(){
    if(!are_OtherPV_loaded_){
      tree_->SetBranchStatus("OtherPV_z", 1); tree_->SetBranchAddress("OtherPV_z", &OtherPV_z_);
      are_OtherPV_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadSv(){
    if(!are_SV_loaded_){
      tree_->SetBranchStatus("SV_dlen", 1); tree_->SetBranchAddress("SV_dlen", &SV_dlen_);
      tree_->SetBranchStatus("SV_dlenSig", 1); tree_->SetBranchAddress("SV_dlenSig", &SV_dlenSig_);
      tree_->SetBranchStatus("SV_pAngle", 1); tree_->SetBranchAddress("SV_pAngle", &SV_pAngle_);
      tree_->SetBranchStatus("SV_chi2", 1); tree_->SetBranchAddress("SV_chi2", &SV_chi2_);
      tree_->SetBranchStatus("SV_eta", 1); tree_->SetBranchAddress("SV_eta", &SV_eta_);
      tree_->SetBranchStatus("SV_mass", 1); tree_->SetBranchAddress("SV_mass", &SV_mass_);
      tree_->SetBranchStatus("SV_ndof", 1); tree_->SetBranchAddress("SV_ndof", &SV_ndof_);
      tree_->SetBranchStatus("SV_phi", 1); tree_->SetBranchAddress("SV_phi", &SV_phi_);
      tree_->SetBranchStatus("SV_pt", 1); tree_->SetBranchAddress("SV_pt", &SV_pt_);
      tree_->SetBranchStatus("SV_x", 1); tree_->SetBranchAddress("SV_x", &SV_x_);
      tree_->SetBranchStatus("SV_y", 1); tree_->SetBranchAddress("SV_y", &SV_y_);
      tree_->SetBranchStatus("SV_z", 1); tree_->SetBranchAddress("SV_z", &SV_z_);
      are_SV_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  

  const Electron Electron(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadElectron();
  
    Electron obj;
    obj.setdeltaEtaSC(Electron_deltaEtaSC_);
    obj.setdr03EcalRecHitSumEt(Electron_dr03EcalRecHitSumEt_);
    obj.setdr03HcalDepth1TowerSumEt(Electron_dr03HcalDepth1TowerSumEt_);
    obj.setdr03TkSumPt(Electron_dr03TkSumPt_);
    obj.setdxy(Electron_dxy_);
    obj.setdxyErr(Electron_dxyErr_);
    obj.setdz(Electron_dz_);
    obj.setdzErr(Electron_dzErr_);
    obj.seteCorr(Electron_eCorr_);
    obj.seteInvMinusPInv(Electron_eInvMinusPInv_);
    obj.setenergyErr(Electron_energyErr_);
    obj.sethoe(Electron_hoe_);
    obj.setip3d(Electron_ip3d_);
    obj.setminiPFRelIso(Electron_miniPFRelIso_all_);
    obj.setminiPFRelIso(Electron_miniPFRelIso_chg_);
    obj.setmvaSpring16GP(Electron_mvaSpring16GP_);
    obj.setmvaSpring16HZZ(Electron_mvaSpring16HZZ_);
    obj.setpfRelIso03(Electron_pfRelIso03_all_);
    obj.setpfRelIso03(Electron_pfRelIso03_chg_);
    obj.setr9(Electron_r9_);
    obj.setsieie(Electron_sieie_);
    obj.setsip3d(Electron_sip3d_);
    obj.setmvaTTH(Electron_mvaTTH_);
    obj.setcharge(Electron_charge_);
    obj.setcutBased(Electron_cutBased_);
    obj.setcutBased(Electron_cutBased_HLTPreSel_);
    obj.setjetIdx(Electron_jetIdx_);
    obj.setpdgId(Electron_pdgId_);
    obj.setphotonIdx(Electron_photonIdx_);
    obj.settightCharge(Electron_tightCharge_);
    obj.setvidNestedWPBitmap(Electron_vidNestedWPBitmap_);
    obj.setconvVeto(Electron_convVeto_);
    obj.setcutBased(Electron_cutBased_HEEP_);
    obj.setisPFcand(Electron_isPFcand_);
    obj.setlostHits(Electron_lostHits_);
    obj.setmvaSpring16GP(Electron_mvaSpring16GP_WP80_);
    obj.setmvaSpring16GP(Electron_mvaSpring16GP_WP90_);
    obj.setmvaSpring16HZZ(Electron_mvaSpring16HZZ_WPL_);
    obj.setgenPartIdx(Electron_genPartIdx_);
    obj.setgenPartFlav(Electron_genPartFlav_);
    obj.setcleanmask(Electron_cleanmask_);
  
    return obj;
  }
  
  const Fatjet FatJet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadFatjet();
  
    Fatjet obj;
    obj.setarea(FatJet_area_);
    obj.setbtagCMVA(FatJet_btagCMVA_);
    obj.setbtagCSVV2(FatJet_btagCSVV2_);
    obj.setbtagDeepB(FatJet_btagDeepB_);
    obj.setbtagHbb(FatJet_btagHbb_);
    obj.setmsoftdrop(FatJet_msoftdrop_);
    obj.setmsoftdrop(FatJet_msoftdrop_chs_);
    obj.setn2b1(FatJet_n2b1_);
    obj.setn3b1(FatJet_n3b1_);
    obj.settau1(FatJet_tau1_);
    obj.settau2(FatJet_tau2_);
    obj.settau3(FatJet_tau3_);
    obj.settau4(FatJet_tau4_);
    obj.setjetId(FatJet_jetId_);
    obj.setsubJetIdx1(FatJet_subJetIdx1_);
    obj.setsubJetIdx2(FatJet_subJetIdx2_);
  
    return obj;
  }
  
  const Genjetak8 GenJetAK8(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGenjetak8();
  
    Genjetak8 obj;
    obj.setpartonFlavour(GenJetAK8_partonFlavour_);
    obj.sethadronFlavour(GenJetAK8_hadronFlavour_);
  
    return obj;
  }
  
  const Genjet GenJet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGenjet();
  
    Genjet obj;
    obj.setpartonFlavour(GenJetAK8_partonFlavour_);
    obj.sethadronFlavour(GenJetAK8_hadronFlavour_);
    obj.setpartonFlavour(GenJet_partonFlavour_);
    obj.sethadronFlavour(GenJet_hadronFlavour_);
  
    return obj;
  }
  
  const Genpart GenPart(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGenpart();
  
    Genpart obj;
    obj.setgenPartIdxMother(GenPart_genPartIdxMother_);
    obj.setpdgId(GenPart_pdgId_);
    obj.setstatus(GenPart_status_);
    obj.setstatusFlags(GenPart_statusFlags_);
  
    return obj;
  }
  
  const Genvistau GenVisTau(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGenvistau();
  
    Genvistau obj;
    obj.setcharge(GenVisTau_charge_);
    obj.setgenPartIdxMother(GenVisTau_genPartIdxMother_);
    obj.setstatus(GenVisTau_status_);
  
    return obj;
  }
  
  const Lhepdfweight LHEPdfWeight(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadLhepdfweight();
  
    Lhepdfweight obj;
    obj.setLHEPdfWeight(LHEPdfWeight);
  
    return obj;
  }
  
  const Lhescaleweight LHEScaleWeight(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadLhescaleweight();
  
    Lhescaleweight obj;
    obj.setLHEScaleWeight(LHEScaleWeight);
  
    return obj;
  }
  
  const Jet Jet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadJet();
  
    Jet obj;
    obj.setarea(Jet_area_);
    obj.setbtagCMVA(Jet_btagCMVA_);
    obj.setbtagCSVV2(Jet_btagCSVV2_);
    obj.setbtagDeepB(Jet_btagDeepB_);
    obj.setbtagDeepC(Jet_btagDeepC_);
    obj.setchEmEF(Jet_chEmEF_);
    obj.setchHEF(Jet_chHEF_);
    obj.setneEmEF(Jet_neEmEF_);
    obj.setneHEF(Jet_neHEF_);
    obj.setqgl(Jet_qgl_);
    obj.setrawFactor(Jet_rawFactor_);
    obj.setbReg(Jet_bReg_);
    obj.setelectronIdx1(Jet_electronIdx1_);
    obj.setelectronIdx2(Jet_electronIdx2_);
    obj.setjetId(Jet_jetId_);
    obj.setmuonIdx1(Jet_muonIdx1_);
    obj.setmuonIdx2(Jet_muonIdx2_);
    obj.setnConstituents(Jet_nConstituents_);
    obj.setnElectrons(Jet_nElectrons_);
    obj.setnMuons(Jet_nMuons_);
    obj.setpuId(Jet_puId_);
    obj.setgenJetIdx(Jet_genJetIdx_);
    obj.sethadronFlavour(Jet_hadronFlavour_);
    obj.setpartonFlavour(Jet_partonFlavour_);
    obj.setcleanmask(Jet_cleanmask_);
  
    return obj;
  }
  
  const Muon Muon(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadMuon();
  
    Muon obj;
    obj.setdxy(Muon_dxy_);
    obj.setdxyErr(Muon_dxyErr_);
    obj.setdz(Muon_dz_);
    obj.setdzErr(Muon_dzErr_);
    obj.setip3d(Muon_ip3d_);
    obj.setminiPFRelIso(Muon_miniPFRelIso_all_);
    obj.setminiPFRelIso(Muon_miniPFRelIso_chg_);
    obj.setpfRelIso03(Muon_pfRelIso03_all_);
    obj.setpfRelIso03(Muon_pfRelIso03_chg_);
    obj.setpfRelIso04(Muon_pfRelIso04_all_);
    obj.setptErr(Muon_ptErr_);
    obj.setsegmentComp(Muon_segmentComp_);
    obj.setsip3d(Muon_sip3d_);
    obj.setmvaTTH(Muon_mvaTTH_);
    obj.setcharge(Muon_charge_);
    obj.setjetIdx(Muon_jetIdx_);
    obj.setnStations(Muon_nStations_);
    obj.setnTrackerLayers(Muon_nTrackerLayers_);
    obj.setpdgId(Muon_pdgId_);
    obj.settightCharge(Muon_tightCharge_);
    obj.sethighPtId(Muon_highPtId_);
    obj.setisPFcand(Muon_isPFcand_);
    obj.setmediumId(Muon_mediumId_);
    obj.setsoftId(Muon_softId_);
    obj.settightId(Muon_tightId_);
    obj.setgenPartIdx(Muon_genPartIdx_);
    obj.setgenPartFlav(Muon_genPartFlav_);
    obj.setcleanmask(Muon_cleanmask_);
  
    return obj;
  }
  
  const Photon Photon(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadPhoton();
  
    Photon obj;
    obj.seteCorr(Photon_eCorr_);
    obj.setenergyErr(Photon_energyErr_);
    obj.sethoe(Photon_hoe_);
    obj.setmvaID(Photon_mvaID_);
    obj.setpfRelIso03(Photon_pfRelIso03_all_);
    obj.setpfRelIso03(Photon_pfRelIso03_chg_);
    obj.setr9(Photon_r9_);
    obj.setsieie(Photon_sieie_);
    obj.setcharge(Photon_charge_);
    obj.setcutBased(Photon_cutBased_);
    obj.setelectronIdx(Photon_electronIdx_);
    obj.setjetIdx(Photon_jetIdx_);
    obj.setpdgId(Photon_pdgId_);
    obj.setvidNestedWPBitmap(Photon_vidNestedWPBitmap_);
    obj.setelectronVeto(Photon_electronVeto_);
    obj.setmvaID(Photon_mvaID_WP80_);
    obj.setmvaID(Photon_mvaID_WP90_);
    obj.setpixelSeed(Photon_pixelSeed_);
    obj.setgenPartIdx(Photon_genPartIdx_);
    obj.setgenPartFlav(Photon_genPartFlav_);
    obj.setcleanmask(Photon_cleanmask_);
  
    return obj;
  }
  
  const Gendressedlepton GenDressedLepton(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGendressedlepton();
  
    Gendressedlepton obj;
    obj.setpdgId(GenDressedLepton_pdgId_);
  
    return obj;
  }
  
  const Softactivityjet SoftActivityJet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadSoftactivityjet();
  
    Softactivityjet obj;
    obj.setSoftActivityJetHT(SoftActivityJetHT);
    obj.setSoftActivityJetHT10(SoftActivityJetHT10);
    obj.setSoftActivityJetHT2(SoftActivityJetHT2);
    obj.setSoftActivityJetHT5(SoftActivityJetHT5);
    obj.setSoftActivityJetNjets10(SoftActivityJetNjets10);
    obj.setSoftActivityJetNjets2(SoftActivityJetNjets2);
    obj.setSoftActivityJetNjets5(SoftActivityJetNjets5);
  
    return obj;
  }
  
  const Subjet SubJet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadSubjet();
  
    Subjet obj;
    obj.setbtagCMVA(SubJet_btagCMVA_);
    obj.setbtagCSVV2(SubJet_btagCSVV2_);
    obj.setbtagDeepB(SubJet_btagDeepB_);
    obj.setn2b1(SubJet_n2b1_);
    obj.setn3b1(SubJet_n3b1_);
    obj.settau1(SubJet_tau1_);
    obj.settau2(SubJet_tau2_);
    obj.settau3(SubJet_tau3_);
    obj.settau4(SubJet_tau4_);
  
    return obj;
  }
  
  const Tau Tau(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadTau();
  
    Tau obj;
    obj.setchargedIso(Tau_chargedIso_);
    obj.setdxy(Tau_dxy_);
    obj.setdz(Tau_dz_);
    obj.setfootprintCorr(Tau_footprintCorr_);
    obj.setleadTkDeltaEta(Tau_leadTkDeltaEta_);
    obj.setleadTkDeltaPhi(Tau_leadTkDeltaPhi_);
    obj.setleadTkPtOverTauPt(Tau_leadTkPtOverTauPt_);
    obj.setneutralIso(Tau_neutralIso_);
    obj.setphotonsOutsideSignalCone(Tau_photonsOutsideSignalCone_);
    obj.setpuCorr(Tau_puCorr_);
    obj.setrawAntiEle(Tau_rawAntiEle_);
    obj.setrawIso(Tau_rawIso_);
    obj.setrawMVAnewDM(Tau_rawMVAnewDM_);
    obj.setrawMVAoldDM(Tau_rawMVAoldDM_);
    obj.setrawMVAoldDMdR03(Tau_rawMVAoldDMdR03_);
    obj.setcharge(Tau_charge_);
    obj.setdecayMode(Tau_decayMode_);
    obj.setjetIdx(Tau_jetIdx_);
    obj.setrawAntiEleCat(Tau_rawAntiEleCat_);
    obj.setidAntiEle(Tau_idAntiEle_);
    obj.setidAntiMu(Tau_idAntiMu_);
    obj.setidDecayMode(Tau_idDecayMode_);
    obj.setidDecayModeNewDMs(Tau_idDecayModeNewDMs_);
    obj.setidMVAnewDM(Tau_idMVAnewDM_);
    obj.setidMVAoldDM(Tau_idMVAoldDM_);
    obj.setidMVAoldDMdR03(Tau_idMVAoldDMdR03_);
    obj.setcleanmask(Tau_cleanmask_);
    obj.setgenPartIdx(Tau_genPartIdx_);
    obj.setgenPartFlav(Tau_genPartFlav_);
  
    return obj;
  }
  
  const Trigobj TrigObj(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadTrigobj();
  
    Trigobj obj;
    obj.setl1pt(TrigObj_l1pt_);
    obj.setl1pt(TrigObj_l1pt_2_);
    obj.setl2pt(TrigObj_l2pt_);
    obj.setid(TrigObj_id_);
    obj.setl1iso(TrigObj_l1iso_);
    obj.setl1charge(TrigObj_l1charge_);
    obj.setfilterBits(TrigObj_filterBits_);
  
    return obj;
  }
  
  const Otherpv OtherPV(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadOtherpv();
  
    Otherpv obj;
    obj.setz(OtherPV_z_);
  
    return obj;
  }
  
  const Sv SV(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadSv();
  
    Sv obj;
    obj.setdlen(SV_dlen_);
    obj.setdlenSig(SV_dlenSig_);
    obj.setpAngle(SV_pAngle_);
    obj.setchi2(SV_chi2_);
    obj.setndof(SV_ndof_);
    obj.setx(SV_x_);
    obj.sety(SV_y_);
    obj.setz(SV_z_);
  
    return obj;
  }
  

private:
  TTree *tree_;
  Long64_t entries_;
  Long64_t current_entry_;
  Float_t CaloMET_phi_;
  Float_t CaloMET_pt_;
  Float_t CaloMET_sumEt_;
  Float_t Electron_deltaEtaSC_;
  Float_t Electron_dr03EcalRecHitSumEt_;
  Float_t Electron_dr03HcalDepth1TowerSumEt_;
  Float_t Electron_dr03TkSumPt_;
  Float_t Electron_dxy_;
  Float_t Electron_dxyErr_;
  Float_t Electron_dz_;
  Float_t Electron_dzErr_;
  Float_t Electron_eCorr_;
  Float_t Electron_eInvMinusPInv_;
  Float_t Electron_energyErr_;
  Float_t Electron_eta_;
  Float_t Electron_hoe_;
  Float_t Electron_ip3d_;
  Float_t Electron_mass_;
  Float_t Electron_miniPFRelIso_all_;
  Float_t Electron_miniPFRelIso_chg_;
  Float_t Electron_mvaSpring16GP_;
  Float_t Electron_mvaSpring16HZZ_;
  Float_t Electron_pfRelIso03_all_;
  Float_t Electron_pfRelIso03_chg_;
  Float_t Electron_phi_;
  Float_t Electron_pt_;
  Float_t Electron_r9_;
  Float_t Electron_sieie_;
  Float_t Electron_sip3d_;
  Float_t Electron_mvaTTH_;
  Int_t Electron_charge_;
  Int_t Electron_cutBased_;
  Int_t Electron_cutBased_HLTPreSel_;
  Int_t Electron_jetIdx_;
  Int_t Electron_pdgId_;
  Int_t Electron_photonIdx_;
  Int_t Electron_tightCharge_;
  Int_t Electron_vidNestedWPBitmap_;
  Bool_t Electron_convVeto_;
  Bool_t Electron_cutBased_HEEP_;
  Bool_t Electron_isPFcand_;
  UChar_t Electron_lostHits_;
  Bool_t Electron_mvaSpring16GP_WP80_;
  Bool_t Electron_mvaSpring16GP_WP90_;
  Bool_t Electron_mvaSpring16HZZ_WPL_;
  Float_t FatJet_area_;
  Float_t FatJet_btagCMVA_;
  Float_t FatJet_btagCSVV2_;
  Float_t FatJet_btagDeepB_;
  Float_t FatJet_btagHbb_;
  Float_t FatJet_eta_;
  Float_t FatJet_mass_;
  Float_t FatJet_msoftdrop_;
  Float_t FatJet_msoftdrop_chs_;
  Float_t FatJet_n2b1_;
  Float_t FatJet_n3b1_;
  Float_t FatJet_phi_;
  Float_t FatJet_pt_;
  Float_t FatJet_tau1_;
  Float_t FatJet_tau2_;
  Float_t FatJet_tau3_;
  Float_t FatJet_tau4_;
  Int_t FatJet_jetId_;
  Int_t FatJet_subJetIdx1_;
  Int_t FatJet_subJetIdx2_;
  Float_t GenJetAK8_eta_;
  Float_t GenJetAK8_mass_;
  Float_t GenJetAK8_phi_;
  Float_t GenJetAK8_pt_;
  Float_t GenJet_eta_;
  Float_t GenJet_mass_;
  Float_t GenJet_phi_;
  Float_t GenJet_pt_;
  Float_t GenPart_eta_;
  Float_t GenPart_mass_;
  Float_t GenPart_phi_;
  Float_t GenPart_pt_;
  Int_t GenPart_genPartIdxMother_;
  Int_t GenPart_pdgId_;
  Int_t GenPart_status_;
  Int_t GenPart_statusFlags_;
  Float_t Generator_binvar_;
  Float_t Generator_scalePDF_;
  Float_t Generator_weight_;
  Float_t Generator_x1_;
  Float_t Generator_x2_;
  Float_t Generator_xpdf1_;
  Float_t Generator_xpdf2_;
  Int_t Generator_id1_;
  Int_t Generator_id2_;
  Float_t GenVisTau_eta_;
  Float_t GenVisTau_mass_;
  Float_t GenVisTau_phi_;
  Float_t GenVisTau_pt_;
  Int_t GenVisTau_charge_;
  Int_t GenVisTau_genPartIdxMother_;
  Int_t GenVisTau_status_;
  Float_t LHEWeight_originalXWGTUP_;
  Float_t Jet_area_;
  Float_t Jet_btagCMVA_;
  Float_t Jet_btagCSVV2_;
  Float_t Jet_btagDeepB_;
  Float_t Jet_btagDeepC_;
  Float_t Jet_chEmEF_;
  Float_t Jet_chHEF_;
  Float_t Jet_eta_;
  Float_t Jet_mass_;
  Float_t Jet_neEmEF_;
  Float_t Jet_neHEF_;
  Float_t Jet_phi_;
  Float_t Jet_pt_;
  Float_t Jet_qgl_;
  Float_t Jet_rawFactor_;
  Float_t Jet_bReg_;
  Int_t Jet_electronIdx1_;
  Int_t Jet_electronIdx2_;
  Int_t Jet_jetId_;
  Int_t Jet_muonIdx1_;
  Int_t Jet_muonIdx2_;
  Int_t Jet_nConstituents_;
  Int_t Jet_nElectrons_;
  Int_t Jet_nMuons_;
  Int_t Jet_puId_;
  Float_t LHE_HT_;
  Float_t LHE_HTIncoming_;
  Float_t LHE_Vpt_;
  UChar_t LHE_Njets_;
  UChar_t LHE_Nb_;
  UChar_t LHE_Nc_;
  UChar_t LHE_Nuds_;
  UChar_t LHE_Nglu_;
  UChar_t LHE_NpNLO_;
  UChar_t LHE_NpLO_;
  Float_t GenMET_phi_;
  Float_t GenMET_pt_;
  Float_t MET_MetUnclustEnUpDeltaX_;
  Float_t MET_MetUnclustEnUpDeltaY_;
  Float_t MET_covXX_;
  Float_t MET_covXY_;
  Float_t MET_covYY_;
  Float_t MET_phi_;
  Float_t MET_pt_;
  Float_t MET_significance_;
  Float_t MET_sumEt_;
  Float_t Muon_dxy_;
  Float_t Muon_dxyErr_;
  Float_t Muon_dz_;
  Float_t Muon_dzErr_;
  Float_t Muon_eta_;
  Float_t Muon_ip3d_;
  Float_t Muon_mass_;
  Float_t Muon_miniPFRelIso_all_;
  Float_t Muon_miniPFRelIso_chg_;
  Float_t Muon_pfRelIso03_all_;
  Float_t Muon_pfRelIso03_chg_;
  Float_t Muon_pfRelIso04_all_;
  Float_t Muon_phi_;
  Float_t Muon_pt_;
  Float_t Muon_ptErr_;
  Float_t Muon_segmentComp_;
  Float_t Muon_sip3d_;
  Float_t Muon_mvaTTH_;
  Int_t Muon_charge_;
  Int_t Muon_jetIdx_;
  Int_t Muon_nStations_;
  Int_t Muon_nTrackerLayers_;
  Int_t Muon_pdgId_;
  Int_t Muon_tightCharge_;
  UChar_t Muon_highPtId_;
  Bool_t Muon_isPFcand_;
  Bool_t Muon_mediumId_;
  Bool_t Muon_softId_;
  Bool_t Muon_tightId_;
  Float_t Photon_eCorr_;
  Float_t Photon_energyErr_;
  Float_t Photon_eta_;
  Float_t Photon_hoe_;
  Float_t Photon_mass_;
  Float_t Photon_mvaID_;
  Float_t Photon_pfRelIso03_all_;
  Float_t Photon_pfRelIso03_chg_;
  Float_t Photon_phi_;
  Float_t Photon_pt_;
  Float_t Photon_r9_;
  Float_t Photon_sieie_;
  Int_t Photon_charge_;
  Int_t Photon_cutBased_;
  Int_t Photon_electronIdx_;
  Int_t Photon_jetIdx_;
  Int_t Photon_pdgId_;
  Int_t Photon_vidNestedWPBitmap_;
  Bool_t Photon_electronVeto_;
  Bool_t Photon_mvaID_WP80_;
  Bool_t Photon_mvaID_WP90_;
  Bool_t Photon_pixelSeed_;
  Float_t Pileup_nTrueInt_;
  Int_t Pileup_nPU_;
  Int_t Pileup_sumEOOT_;
  Int_t Pileup_sumLOOT_;
  Float_t PuppiMET_phi_;
  Float_t PuppiMET_pt_;
  Float_t PuppiMET_sumEt_;
  Float_t RawMET_phi_;
  Float_t RawMET_pt_;
  Float_t RawMET_sumEt_;
  Float_t GenDressedLepton_eta_;
  Float_t GenDressedLepton_mass_;
  Float_t GenDressedLepton_phi_;
  Float_t GenDressedLepton_pt_;
  Int_t GenDressedLepton_pdgId_;
  Float_t SoftActivityJet_eta_;
  Float_t SoftActivityJet_phi_;
  Float_t SoftActivityJet_pt_;
  Float_t SubJet_btagCMVA_;
  Float_t SubJet_btagCSVV2_;
  Float_t SubJet_btagDeepB_;
  Float_t SubJet_eta_;
  Float_t SubJet_mass_;
  Float_t SubJet_n2b1_;
  Float_t SubJet_n3b1_;
  Float_t SubJet_phi_;
  Float_t SubJet_pt_;
  Float_t SubJet_tau1_;
  Float_t SubJet_tau2_;
  Float_t SubJet_tau3_;
  Float_t SubJet_tau4_;
  Float_t Tau_chargedIso_;
  Float_t Tau_dxy_;
  Float_t Tau_dz_;
  Float_t Tau_eta_;
  Float_t Tau_footprintCorr_;
  Float_t Tau_leadTkDeltaEta_;
  Float_t Tau_leadTkDeltaPhi_;
  Float_t Tau_leadTkPtOverTauPt_;
  Float_t Tau_mass_;
  Float_t Tau_neutralIso_;
  Float_t Tau_phi_;
  Float_t Tau_photonsOutsideSignalCone_;
  Float_t Tau_pt_;
  Float_t Tau_puCorr_;
  Float_t Tau_rawAntiEle_;
  Float_t Tau_rawIso_;
  Float_t Tau_rawMVAnewDM_;
  Float_t Tau_rawMVAoldDM_;
  Float_t Tau_rawMVAoldDMdR03_;
  Int_t Tau_charge_;
  Int_t Tau_decayMode_;
  Int_t Tau_jetIdx_;
  Int_t Tau_rawAntiEleCat_;
  UChar_t Tau_idAntiEle_;
  UChar_t Tau_idAntiMu_;
  Bool_t Tau_idDecayMode_;
  Bool_t Tau_idDecayModeNewDMs_;
  UChar_t Tau_idMVAnewDM_;
  UChar_t Tau_idMVAoldDM_;
  UChar_t Tau_idMVAoldDMdR03_;
  Float_t TkMET_phi_;
  Float_t TkMET_pt_;
  Float_t TkMET_sumEt_;
  Float_t TrigObj_pt_;
  Float_t TrigObj_eta_;
  Float_t TrigObj_phi_;
  Float_t TrigObj_l1pt_;
  Float_t TrigObj_l1pt_2_;
  Float_t TrigObj_l2pt_;
  Int_t TrigObj_id_;
  Int_t TrigObj_l1iso_;
  Int_t TrigObj_l1charge_;
  Int_t TrigObj_filterBits_;
  Float_t OtherPV_z_;
  Float_t PV_ndof_;
  Float_t PV_x_;
  Float_t PV_y_;
  Float_t PV_z_;
  Float_t PV_chi2_;
  Float_t PV_score_;
  Int_t PV_npvs_;
  Int_t PV_npvsGood_;
  Float_t SV_dlen_;
  Float_t SV_dlenSig_;
  Float_t SV_pAngle_;
  Int_t Electron_genPartIdx_;
  UChar_t Electron_genPartFlav_;
  Int_t GenJetAK8_partonFlavour_;
  UChar_t GenJetAK8_hadronFlavour_;
  Int_t GenJet_partonFlavour_;
  UChar_t GenJet_hadronFlavour_;
  Int_t Jet_genJetIdx_;
  Int_t Jet_hadronFlavour_;
  Int_t Jet_partonFlavour_;
  Int_t Muon_genPartIdx_;
  UChar_t Muon_genPartFlav_;
  Int_t Photon_genPartIdx_;
  UChar_t Photon_genPartFlav_;
  Float_t MET_fiducialGenPhi_;
  Float_t MET_fiducialGenPt_;
  UChar_t Electron_cleanmask_;
  UChar_t Jet_cleanmask_;
  UChar_t Muon_cleanmask_;
  UChar_t Photon_cleanmask_;
  UChar_t Tau_cleanmask_;
  Float_t SV_chi2_;
  Float_t SV_eta_;
  Float_t SV_mass_;
  Float_t SV_ndof_;
  Float_t SV_phi_;
  Float_t SV_pt_;
  Float_t SV_x_;
  Float_t SV_y_;
  Float_t SV_z_;
  Int_t Tau_genPartIdx_;
  UChar_t Tau_genPartFlav_;
  Bool_t L1simulation_step_;
  Bool_t HLT_AK8PFJet360_TrimMass30_;
  Bool_t HLT_AK8PFJet400_TrimMass30_;
  Bool_t HLT_AK8PFHT750_TrimMass50_;
  Bool_t HLT_AK8PFHT800_TrimMass50_;
  Bool_t HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_;
  Bool_t HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_;
  Bool_t HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_;
  Bool_t HLT_AK8DiPFJet300_200_TrimMass30_;
  Bool_t HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_;
  Bool_t HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_;
  Bool_t HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_;
  Bool_t HLT_AK8DiPFJet280_200_TrimMass30_;
  Bool_t HLT_AK8DiPFJet250_200_TrimMass30_;
  Bool_t HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_;
  Bool_t HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_;
  Bool_t HLT_CaloJet260_;
  Bool_t HLT_CaloJet500_NoJetID_;
  Bool_t HLT_Dimuon13_PsiPrime_;
  Bool_t HLT_Dimuon13_Upsilon_;
  Bool_t HLT_Dimuon20_Jpsi_;
  Bool_t HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_;
  Bool_t HLT_DoubleEle25_CaloIdL_GsfTrkIdVL_;
  Bool_t HLT_DoubleEle33_CaloIdL_;
  Bool_t HLT_DoubleEle33_CaloIdL_MW_;
  Bool_t HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_;
  Bool_t HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_;
  Bool_t HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_;
  Bool_t HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_;
  Bool_t HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_;
  Bool_t HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_;
  Bool_t HLT_DoubleMu33NoFiltersNoVtx_;
  Bool_t HLT_DoubleMu38NoFiltersNoVtx_;
  Bool_t HLT_DoubleMu23NoFiltersNoVtxDisplaced_;
  Bool_t HLT_DoubleMu28NoFiltersNoVtxDisplaced_;
  Bool_t HLT_DoubleMu0_;
  Bool_t HLT_DoubleMu4_3_Bs_;
  Bool_t HLT_DoubleMu4_3_Jpsi_Displaced_;
  Bool_t HLT_DoubleMu4_JpsiTrk_Displaced_;
  Bool_t HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_;
  Bool_t HLT_DoubleMu3_Trk_Tau3mu_;
  Bool_t HLT_DoubleMu4_PsiPrimeTrk_Displaced_;
  Bool_t HLT_Mu7p5_L2Mu2_Jpsi_;
  Bool_t HLT_Mu7p5_L2Mu2_Upsilon_;
  Bool_t HLT_Mu7p5_Track2_Jpsi_;
  Bool_t HLT_Mu7p5_Track3p5_Jpsi_;
  Bool_t HLT_Mu7p5_Track7_Jpsi_;
  Bool_t HLT_Mu7p5_Track2_Upsilon_;
  Bool_t HLT_Mu7p5_Track3p5_Upsilon_;
  Bool_t HLT_Mu7p5_Track7_Upsilon_;
  Bool_t HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_;
  Bool_t HLT_Dimuon0er16_Jpsi_NoVertexing_;
  Bool_t HLT_Dimuon6_Jpsi_NoVertexing_;
  Bool_t HLT_Photon150_;
  Bool_t HLT_Photon90_CaloIdL_HT300_;
  Bool_t HLT_HT250_CaloMET70_;
  Bool_t HLT_DoublePhoton60_;
  Bool_t HLT_DoublePhoton85_;
  Bool_t HLT_Ele17_Ele8_Gsf_;
  Bool_t HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28_;
  Bool_t HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29_;
  Bool_t HLT_Ele22_eta2p1_WPLoose_Gsf_;
  Bool_t HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_;
  Bool_t HLT_Ele23_WPLoose_Gsf_;
  Bool_t HLT_Ele23_WPLoose_Gsf_WHbbBoost_;
  Bool_t HLT_Ele24_eta2p1_WPLoose_Gsf_;
  Bool_t HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_;
  Bool_t HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_;
  Bool_t HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_;
  Bool_t HLT_Ele25_WPTight_Gsf_;
  Bool_t HLT_Ele25_eta2p1_WPLoose_Gsf_;
  Bool_t HLT_Ele25_eta2p1_WPTight_Gsf_;
  Bool_t HLT_Ele27_WPLoose_Gsf_;
  Bool_t HLT_Ele27_WPLoose_Gsf_WHbbBoost_;
  Bool_t HLT_Ele27_WPTight_Gsf_;
  Bool_t HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_;
  Bool_t HLT_Ele27_eta2p1_WPLoose_Gsf_;
  Bool_t HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_;
  Bool_t HLT_Ele27_eta2p1_WPTight_Gsf_;
  Bool_t HLT_Ele30_WPTight_Gsf_;
  Bool_t HLT_Ele30_eta2p1_WPLoose_Gsf_;
  Bool_t HLT_Ele30_eta2p1_WPTight_Gsf_;
  Bool_t HLT_Ele32_WPTight_Gsf_;
  Bool_t HLT_Ele32_eta2p1_WPLoose_Gsf_;
  Bool_t HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_;
  Bool_t HLT_Ele32_eta2p1_WPTight_Gsf_;
  Bool_t HLT_Ele35_WPLoose_Gsf_;
  Bool_t HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_;
  Bool_t HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_;
  Bool_t HLT_Ele45_WPLoose_Gsf_;
  Bool_t HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded_;
  Bool_t HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_;
  Bool_t HLT_Ele105_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_Ele30WP60_SC4_Mass55_;
  Bool_t HLT_Ele30WP60_Ele8_Mass55_;
  Bool_t HLT_HT200_;
  Bool_t HLT_HT275_;
  Bool_t HLT_HT325_;
  Bool_t HLT_HT425_;
  Bool_t HLT_HT575_;
  Bool_t HLT_HT410to430_;
  Bool_t HLT_HT430to450_;
  Bool_t HLT_HT450to470_;
  Bool_t HLT_HT470to500_;
  Bool_t HLT_HT500to550_;
  Bool_t HLT_HT550to650_;
  Bool_t HLT_HT650_;
  Bool_t HLT_Mu16_eta2p1_MET30_;
  Bool_t HLT_IsoMu16_eta2p1_MET30_;
  Bool_t HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1_;
  Bool_t HLT_IsoMu17_eta2p1_;
  Bool_t HLT_IsoMu17_eta2p1_LooseIsoPFTau20_;
  Bool_t HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_;
  Bool_t HLT_DoubleIsoMu17_eta2p1_;
  Bool_t HLT_DoubleIsoMu17_eta2p1_noDzCut_;
  Bool_t HLT_IsoMu18_;
  Bool_t HLT_IsoMu19_eta2p1_LooseIsoPFTau20_;
  Bool_t HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_;
  Bool_t HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_;
  Bool_t HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_;
  Bool_t HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_;
  Bool_t HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_;
  Bool_t HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_;
  Bool_t HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_;
  Bool_t HLT_IsoMu20_;
  Bool_t HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_;
  Bool_t HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1_;
  Bool_t HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_;
  Bool_t HLT_IsoMu22_;
  Bool_t HLT_IsoMu22_eta2p1_;
  Bool_t HLT_IsoMu24_;
  Bool_t HLT_IsoMu27_;
  Bool_t HLT_IsoTkMu18_;
  Bool_t HLT_IsoTkMu20_;
  Bool_t HLT_IsoTkMu22_;
  Bool_t HLT_IsoTkMu22_eta2p1_;
  Bool_t HLT_IsoTkMu24_;
  Bool_t HLT_IsoTkMu27_;
  Bool_t HLT_JetE30_NoBPTX3BX_;
  Bool_t HLT_JetE30_NoBPTX_;
  Bool_t HLT_JetE50_NoBPTX3BX_;
  Bool_t HLT_JetE70_NoBPTX3BX_;
  Bool_t HLT_L1SingleMu18_;
  Bool_t HLT_L2Mu10_;
  Bool_t HLT_L1SingleMuOpen_;
  Bool_t HLT_L1SingleMuOpen_DT_;
  Bool_t HLT_L2DoubleMu23_NoVertex_;
  Bool_t HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_;
  Bool_t HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10_;
  Bool_t HLT_L2Mu10_NoVertex_NoBPTX3BX_;
  Bool_t HLT_L2Mu10_NoVertex_NoBPTX_;
  Bool_t HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_;
  Bool_t HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_;
  Bool_t HLT_LooseIsoPFTau50_Trk30_eta2p1_;
  Bool_t HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_;
  Bool_t HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_;
  Bool_t HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_;
  Bool_t HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_;
  Bool_t HLT_PFTau120_eta2p1_;
  Bool_t HLT_PFTau140_eta2p1_;
  Bool_t HLT_VLooseIsoPFTau120_Trk50_eta2p1_;
  Bool_t HLT_VLooseIsoPFTau140_Trk50_eta2p1_;
  Bool_t HLT_Mu17_Mu8_;
  Bool_t HLT_Mu17_Mu8_DZ_;
  Bool_t HLT_Mu17_Mu8_SameSign_;
  Bool_t HLT_Mu17_Mu8_SameSign_DZ_;
  Bool_t HLT_Mu20_Mu10_;
  Bool_t HLT_Mu20_Mu10_DZ_;
  Bool_t HLT_Mu20_Mu10_SameSign_;
  Bool_t HLT_Mu20_Mu10_SameSign_DZ_;
  Bool_t HLT_Mu17_TkMu8_DZ_;
  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_;
  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_;
  Bool_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_;
  Bool_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_;
  Bool_t HLT_Mu25_TkMu0_dEta18_Onia_;
  Bool_t HLT_Mu27_TkMu8_;
  Bool_t HLT_Mu30_TkMu11_;
  Bool_t HLT_Mu30_eta2p1_PFJet150_PFJet50_;
  Bool_t HLT_Mu40_TkMu11_;
  Bool_t HLT_Mu40_eta2p1_PFJet200_PFJet50_;
  Bool_t HLT_Mu20_;
  Bool_t HLT_TkMu17_;
  Bool_t HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_;
  Bool_t HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_;
  Bool_t HLT_TkMu20_;
  Bool_t HLT_Mu24_eta2p1_;
  Bool_t HLT_TkMu24_eta2p1_;
  Bool_t HLT_Mu27_;
  Bool_t HLT_TkMu27_;
  Bool_t HLT_Mu45_eta2p1_;
  Bool_t HLT_Mu50_;
  Bool_t HLT_TkMu50_;
  Bool_t HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_;
  Bool_t HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL_;
  Bool_t HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL_;
  Bool_t HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL_;
  Bool_t HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL_;
  Bool_t HLT_DoubleMu18NoFiltersNoVtx_;
  Bool_t HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight_;
  Bool_t HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose_;
  Bool_t HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose_;
  Bool_t HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight_;
  Bool_t HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose_;
  Bool_t HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose_;
  Bool_t HLT_Mu28NoFiltersNoVtx_CentralCaloJet40_;
  Bool_t HLT_PFHT300_PFMET100_;
  Bool_t HLT_PFHT300_PFMET110_;
  Bool_t HLT_PFHT550_4JetPt50_;
  Bool_t HLT_PFHT650_4JetPt50_;
  Bool_t HLT_PFHT750_4JetPt50_;
  Bool_t HLT_PFHT750_4JetPt70_;
  Bool_t HLT_PFHT750_4JetPt80_;
  Bool_t HLT_PFHT800_4JetPt50_;
  Bool_t HLT_PFHT850_4JetPt50_;
  Bool_t HLT_PFJet15_NoCaloMatched_;
  Bool_t HLT_PFJet25_NoCaloMatched_;
  Bool_t HLT_DiPFJet15_NoCaloMatched_;
  Bool_t HLT_DiPFJet25_NoCaloMatched_;
  Bool_t HLT_DiPFJet15_FBEta3_NoCaloMatched_;
  Bool_t HLT_DiPFJet25_FBEta3_NoCaloMatched_;
  Bool_t HLT_DiPFJetAve15_HFJEC_;
  Bool_t HLT_DiPFJetAve25_HFJEC_;
  Bool_t HLT_DiPFJetAve35_HFJEC_;
  Bool_t HLT_AK8PFJet40_;
  Bool_t HLT_AK8PFJet60_;
  Bool_t HLT_AK8PFJet80_;
  Bool_t HLT_AK8PFJet140_;
  Bool_t HLT_AK8PFJet200_;
  Bool_t HLT_AK8PFJet260_;
  Bool_t HLT_AK8PFJet320_;
  Bool_t HLT_AK8PFJet400_;
  Bool_t HLT_AK8PFJet450_;
  Bool_t HLT_AK8PFJet500_;
  Bool_t HLT_PFJet40_;
  Bool_t HLT_PFJet60_;
  Bool_t HLT_PFJet80_;
  Bool_t HLT_PFJet140_;
  Bool_t HLT_PFJet200_;
  Bool_t HLT_PFJet260_;
  Bool_t HLT_PFJet320_;
  Bool_t HLT_PFJet400_;
  Bool_t HLT_PFJet450_;
  Bool_t HLT_PFJet500_;
  Bool_t HLT_DiPFJetAve40_;
  Bool_t HLT_DiPFJetAve60_;
  Bool_t HLT_DiPFJetAve80_;
  Bool_t HLT_DiPFJetAve140_;
  Bool_t HLT_DiPFJetAve200_;
  Bool_t HLT_DiPFJetAve260_;
  Bool_t HLT_DiPFJetAve320_;
  Bool_t HLT_DiPFJetAve400_;
  Bool_t HLT_DiPFJetAve500_;
  Bool_t HLT_DiPFJetAve60_HFJEC_;
  Bool_t HLT_DiPFJetAve80_HFJEC_;
  Bool_t HLT_DiPFJetAve100_HFJEC_;
  Bool_t HLT_DiPFJetAve160_HFJEC_;
  Bool_t HLT_DiPFJetAve220_HFJEC_;
  Bool_t HLT_DiPFJetAve300_HFJEC_;
  Bool_t HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140_;
  Bool_t HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80_;
  Bool_t HLT_DiCentralPFJet170_;
  Bool_t HLT_SingleCentralPFJet170_CFMax0p1_;
  Bool_t HLT_DiCentralPFJet170_CFMax0p1_;
  Bool_t HLT_DiCentralPFJet220_CFMax0p3_;
  Bool_t HLT_DiCentralPFJet330_CFMax0p5_;
  Bool_t HLT_DiCentralPFJet430_;
  Bool_t HLT_PFHT125_;
  Bool_t HLT_PFHT200_;
  Bool_t HLT_PFHT250_;
  Bool_t HLT_PFHT300_;
  Bool_t HLT_PFHT350_;
  Bool_t HLT_PFHT400_;
  Bool_t HLT_PFHT475_;
  Bool_t HLT_PFHT600_;
  Bool_t HLT_PFHT650_;
  Bool_t HLT_PFHT800_;
  Bool_t HLT_PFHT900_;
  Bool_t HLT_PFHT200_PFAlphaT0p51_;
  Bool_t HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57_;
  Bool_t HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_;
  Bool_t HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55_;
  Bool_t HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_;
  Bool_t HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53_;
  Bool_t HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_;
  Bool_t HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52_;
  Bool_t HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53_;
  Bool_t HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51_;
  Bool_t HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52_;
  Bool_t HLT_MET60_IsoTrk35_Loose_;
  Bool_t HLT_MET75_IsoTrk50_;
  Bool_t HLT_MET90_IsoTrk50_;
  Bool_t HLT_PFMET120_BTagCSV_p067_;
  Bool_t HLT_PFMET120_Mu5_;
  Bool_t HLT_PFMET170_NotCleaned_;
  Bool_t HLT_PFMET170_NoiseCleaned_;
  Bool_t HLT_PFMET170_HBHECleaned_;
  Bool_t HLT_PFMET170_JetIdCleaned_;
  Bool_t HLT_PFMET170_BeamHaloCleaned_;
  Bool_t HLT_PFMET170_HBHE_BeamHaloCleaned_;
  Bool_t HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned_;
  Bool_t HLT_PFMET90_PFMHT90_IDTight_;
  Bool_t HLT_PFMET100_PFMHT100_IDTight_;
  Bool_t HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned_;
  Bool_t HLT_PFMET110_PFMHT110_IDTight_;
  Bool_t HLT_PFMET120_PFMHT120_IDTight_;
  Bool_t HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_;
  Bool_t HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_;
  Bool_t HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200_;
  Bool_t HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460_;
  Bool_t HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_;
  Bool_t HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_;
  Bool_t HLT_QuadPFJet_VBF_;
  Bool_t HLT_L1_TripleJet_VBF_;
  Bool_t HLT_QuadJet45_TripleBTagCSV_p087_;
  Bool_t HLT_QuadJet45_DoubleBTagCSV_p087_;
  Bool_t HLT_DoubleJet90_Double30_TripleBTagCSV_p087_;
  Bool_t HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_;
  Bool_t HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_;
  Bool_t HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_;
  Bool_t HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172_;
  Bool_t HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_;
  Bool_t HLT_DoubleJetsC100_SingleBTagCSV_p026_;
  Bool_t HLT_DoubleJetsC100_SingleBTagCSV_p014_;
  Bool_t HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350_;
  Bool_t HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350_;
  Bool_t HLT_Photon135_PFMET100_;
  Bool_t HLT_Photon20_CaloIdVL_IsoL_;
  Bool_t HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40_;
  Bool_t HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_;
  Bool_t HLT_Photon250_NoHE_;
  Bool_t HLT_Photon300_NoHE_;
  Bool_t HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60_;
  Bool_t HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15_;
  Bool_t HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_;
  Bool_t HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_;
  Bool_t HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_;
  Bool_t HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_;
  Bool_t HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40_;
  Bool_t HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_;
  Bool_t HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40_;
  Bool_t HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_;
  Bool_t HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_;
  Bool_t HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_;
  Bool_t HLT_Mu8_TrkIsoVVL_;
  Bool_t HLT_Mu17_TrkIsoVVL_;
  Bool_t HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t HLT_BTagMu_DiJet20_Mu5_;
  Bool_t HLT_BTagMu_DiJet40_Mu5_;
  Bool_t HLT_BTagMu_DiJet70_Mu5_;
  Bool_t HLT_BTagMu_DiJet110_Mu5_;
  Bool_t HLT_BTagMu_DiJet170_Mu5_;
  Bool_t HLT_BTagMu_Jet300_Mu5_;
  Bool_t HLT_BTagMu_AK8Jet300_Mu5_;
  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_;
  Bool_t HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_;
  Bool_t HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_;
  Bool_t HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL_;
  Bool_t HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL_;
  Bool_t HLT_Mu8_DiEle12_CaloIdL_TrackIdL_;
  Bool_t HLT_Mu12_Photon25_CaloIdL_;
  Bool_t HLT_Mu12_Photon25_CaloIdL_L1ISO_;
  Bool_t HLT_Mu12_Photon25_CaloIdL_L1OR_;
  Bool_t HLT_Mu17_Photon22_CaloIdL_L1ISO_;
  Bool_t HLT_Mu17_Photon30_CaloIdL_L1ISO_;
  Bool_t HLT_Mu17_Photon35_CaloIdL_L1ISO_;
  Bool_t HLT_DiMu9_Ele9_CaloIdL_TrackIdL_;
  Bool_t HLT_TripleMu_5_3_3_;
  Bool_t HLT_TripleMu_12_10_5_;
  Bool_t HLT_Mu3er_PFHT140_PFMET125_;
  Bool_t HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067_;
  Bool_t HLT_Mu6_PFHT200_PFMET100_;
  Bool_t HLT_Mu14er_PFMET100_;
  Bool_t HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Ele17_CaloIdL_GsfTrkIdVL_;
  Bool_t HLT_Ele17_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Ele23_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_;
  Bool_t HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_;
  Bool_t HLT_Photon22_;
  Bool_t HLT_Photon30_;
  Bool_t HLT_Photon36_;
  Bool_t HLT_Photon50_;
  Bool_t HLT_Photon75_;
  Bool_t HLT_Photon90_;
  Bool_t HLT_Photon120_;
  Bool_t HLT_Photon175_;
  Bool_t HLT_Photon165_HE10_;
  Bool_t HLT_Photon22_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon30_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon36_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon50_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon75_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon90_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon120_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon165_R9Id90_HE10_IsoM_;
  Bool_t HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_;
  Bool_t HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_;
  Bool_t HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_;
  Bool_t HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55_;
  Bool_t HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_;
  Bool_t HLT_Dimuon0_Jpsi_Muon_;
  Bool_t HLT_Dimuon0_Upsilon_Muon_;
  Bool_t HLT_QuadMuon0_Dimuon0_Jpsi_;
  Bool_t HLT_QuadMuon0_Dimuon0_Upsilon_;
  Bool_t HLT_Rsq0p25_Calo_;
  Bool_t HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo_;
  Bool_t HLT_RsqMR240_Rsq0p09_MR200_Calo_;
  Bool_t HLT_Rsq0p25_;
  Bool_t HLT_Rsq0p30_;
  Bool_t HLT_RsqMR240_Rsq0p09_MR200_;
  Bool_t HLT_RsqMR240_Rsq0p09_MR200_4jet_;
  Bool_t HLT_RsqMR270_Rsq0p09_MR200_;
  Bool_t HLT_RsqMR270_Rsq0p09_MR200_4jet_;
  Bool_t HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200_;
  Bool_t HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_;
  Bool_t HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_;
  Bool_t HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_;
  Bool_t HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_;
  Bool_t HLT_HT200_DisplacedDijet40_DisplacedTrack_;
  Bool_t HLT_HT250_DisplacedDijet40_DisplacedTrack_;
  Bool_t HLT_HT350_DisplacedDijet40_DisplacedTrack_;
  Bool_t HLT_HT350_DisplacedDijet80_DisplacedTrack_;
  Bool_t HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack_;
  Bool_t HLT_HT350_DisplacedDijet40_Inclusive_;
  Bool_t HLT_HT400_DisplacedDijet40_Inclusive_;
  Bool_t HLT_HT500_DisplacedDijet40_Inclusive_;
  Bool_t HLT_HT550_DisplacedDijet40_Inclusive_;
  Bool_t HLT_HT550_DisplacedDijet80_Inclusive_;
  Bool_t HLT_HT650_DisplacedDijet80_Inclusive_;
  Bool_t HLT_HT750_DisplacedDijet80_Inclusive_;
  Bool_t HLT_VBF_DisplacedJet40_DisplacedTrack_;
  Bool_t HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5_;
  Bool_t HLT_VBF_DisplacedJet40_TightID_DisplacedTrack_;
  Bool_t HLT_VBF_DisplacedJet40_Hadronic_;
  Bool_t HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack_;
  Bool_t HLT_VBF_DisplacedJet40_TightID_Hadronic_;
  Bool_t HLT_VBF_DisplacedJet40_VTightID_Hadronic_;
  Bool_t HLT_VBF_DisplacedJet40_VVTightID_Hadronic_;
  Bool_t HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_;
  Bool_t HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_;
  Bool_t HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_;
  Bool_t HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_;
  Bool_t HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_;
  Bool_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_;
  Bool_t HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight_;
  Bool_t HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight_;
  Bool_t HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_;
  Bool_t HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_;
  Bool_t HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_;
  Bool_t HLT_Photon90_CaloIdL_PFHT500_;
  Bool_t HLT_DoubleMu8_Mass8_PFHT250_;
  Bool_t HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250_;
  Bool_t HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250_;
  Bool_t HLT_DoubleMu8_Mass8_PFHT300_;
  Bool_t HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_;
  Bool_t HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_;
  Bool_t HLT_Mu10_CentralPFJet30_BTagCSV_p13_;
  Bool_t HLT_DoubleMu3_PFMET50_;
  Bool_t HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_;
  Bool_t HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT350_PFMET50_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT600_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT350_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT400_PFMET50_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT400_;
  Bool_t HLT_Ele50_IsoVVVL_PFHT400_;
  Bool_t HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_;
  Bool_t HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_;
  Bool_t HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT350_PFMET50_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT600_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT350_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT400_PFMET50_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT400_;
  Bool_t HLT_Mu50_IsoVVVL_PFHT400_;
  Bool_t HLT_Dimuon16_Jpsi_;
  Bool_t HLT_Dimuon10_Jpsi_Barrel_;
  Bool_t HLT_Dimuon8_PsiPrime_Barrel_;
  Bool_t HLT_Dimuon8_Upsilon_Barrel_;
  Bool_t HLT_Dimuon0_Phi_Barrel_;
  Bool_t HLT_Mu16_TkMu0_dEta18_Onia_;
  Bool_t HLT_Mu16_TkMu0_dEta18_Phi_;
  Bool_t HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_;
  Bool_t HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_;
  Bool_t HLT_Mu8_;
  Bool_t HLT_Mu17_;
  Bool_t HLT_Mu3_PFJet40_;
  Bool_t HLT_Ele8_CaloIdM_TrackIdM_PFJet30_;
  Bool_t HLT_Ele12_CaloIdM_TrackIdM_PFJet30_;
  Bool_t HLT_Ele17_CaloIdM_TrackIdM_PFJet30_;
  Bool_t HLT_Ele23_CaloIdM_TrackIdM_PFJet30_;
  Bool_t HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_;
  Bool_t HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_;
  Bool_t HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_;
  Bool_t HLT_PFHT450_SixJet40_BTagCSV_p056_;
  Bool_t HLT_PFHT400_SixJet30_;
  Bool_t HLT_PFHT450_SixJet40_;
  Bool_t HLT_Ele115_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_Mu55_;
  Bool_t HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_;
  Bool_t HLT_Photon90_CaloIdL_PFHT600_;
  Bool_t HLT_PixelTracks_Multiplicity60ForEndOfFill_;
  Bool_t HLT_PixelTracks_Multiplicity85ForEndOfFill_;
  Bool_t HLT_PixelTracks_Multiplicity110ForEndOfFill_;
  Bool_t HLT_PixelTracks_Multiplicity135ForEndOfFill_;
  Bool_t HLT_PixelTracks_Multiplicity160ForEndOfFill_;
  Bool_t HLT_FullTracks_Multiplicity80_;
  Bool_t HLT_FullTracks_Multiplicity100_;
  Bool_t HLT_FullTracks_Multiplicity130_;
  Bool_t HLT_FullTracks_Multiplicity150_;
  Bool_t HLT_ECALHT800_;
  Bool_t HLT_DiSC30_18_EIso_AND_HE_Mass70_;
  Bool_t HLT_Photon125_;
  Bool_t HLT_MET100_;
  Bool_t HLT_MET150_;
  Bool_t HLT_MET200_;
  Bool_t HLT_Ele27_HighEta_Ele20_Mass55_;
  Bool_t HLT_L1FatEvents_;
  Bool_t HLT_Physics_;
  Bool_t HLT_L1FatEvents_part0_;
  Bool_t HLT_L1FatEvents_part1_;
  Bool_t HLT_L1FatEvents_part2_;
  Bool_t HLT_L1FatEvents_part3_;
  Bool_t HLT_Random_;
  Bool_t HLT_ZeroBias_;
  Bool_t HLT_AK4CaloJet30_;
  Bool_t HLT_AK4CaloJet40_;
  Bool_t HLT_AK4CaloJet50_;
  Bool_t HLT_AK4CaloJet80_;
  Bool_t HLT_AK4CaloJet100_;
  Bool_t HLT_AK4PFJet30_;
  Bool_t HLT_AK4PFJet50_;
  Bool_t HLT_AK4PFJet80_;
  Bool_t HLT_AK4PFJet100_;
  Bool_t HLT_HISinglePhoton10_;
  Bool_t HLT_HISinglePhoton15_;
  Bool_t HLT_HISinglePhoton20_;
  Bool_t HLT_HISinglePhoton40_;
  Bool_t HLT_HISinglePhoton60_;
  Bool_t HLT_EcalCalibration_;
  Bool_t HLT_HcalCalibration_;
  Bool_t HLT_GlobalRunHPDNoise_;
  Bool_t HLT_L1BptxMinus_;
  Bool_t HLT_L1BptxPlus_;
  Bool_t HLT_L1NotBptxOR_;
  Bool_t HLT_L1BeamGasMinus_;
  Bool_t HLT_L1BeamGasPlus_;
  Bool_t HLT_L1BptxXOR_;
  Bool_t HLT_L1MinimumBiasHF_OR_;
  Bool_t HLT_L1MinimumBiasHF_AND_;
  Bool_t HLT_HcalNZS_;
  Bool_t HLT_HcalPhiSym_;
  Bool_t HLT_HcalIsolatedbunch_;
  Bool_t HLT_ZeroBias_FirstCollisionAfterAbortGap_;
  Bool_t HLT_ZeroBias_FirstCollisionAfterAbortGap_copy_;
  Bool_t HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS_;
  Bool_t HLT_ZeroBias_IsolatedBunches_;
  Bool_t HLT_ZeroBias_FirstCollisionInTrain_;
  Bool_t HLT_ZeroBias_FirstBXAfterTrain_;
  Bool_t HLT_Photon500_;
  Bool_t HLT_Photon600_;
  Bool_t HLT_Mu300_;
  Bool_t HLT_Mu350_;
  Bool_t HLT_MET250_;
  Bool_t HLT_MET300_;
  Bool_t HLT_MET600_;
  Bool_t HLT_MET700_;
  Bool_t HLT_PFMET300_;
  Bool_t HLT_PFMET400_;
  Bool_t HLT_PFMET500_;
  Bool_t HLT_PFMET600_;
  Bool_t HLT_Ele250_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_Ele300_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_HT2000_;
  Bool_t HLT_HT2500_;
  Bool_t HLT_IsoTrackHE_;
  Bool_t HLT_IsoTrackHB_;
  Bool_t Flag_HBHENoiseFilter_;
  Bool_t Flag_HBHENoiseIsoFilter_;
  Bool_t Flag_CSCTightHaloFilter_;
  Bool_t Flag_CSCTightHaloTrkMuUnvetoFilter_;
  Bool_t Flag_CSCTightHalo2015Filter_;
  Bool_t Flag_globalTightHalo2016Filter_;
  Bool_t Flag_globalSuperTightHalo2016Filter_;
  Bool_t Flag_HcalStripHaloFilter_;
  Bool_t Flag_hcalLaserEventFilter_;
  Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter_;
  Bool_t Flag_EcalDeadCellBoundaryEnergyFilter_;
  Bool_t Flag_goodVertices_;
  Bool_t Flag_eeBadScFilter_;
  Bool_t Flag_ecalLaserCorrFilter_;
  Bool_t Flag_trkPOGFilters_;
  Bool_t Flag_chargedHadronTrackResolutionFilter_;
  Bool_t Flag_muonBadTrackFilter_;
  Bool_t Flag_trkPOG_manystripclus53X_;
  Bool_t Flag_trkPOG_toomanystripclus53X_;
  Bool_t Flag_trkPOG_logErrorTooManyClusters_;
  Bool_t Flag_METFilters_;
  bool are_Electron_loaded_;
  Electron Electron_;
  bool are_FatJet_loaded_;
  Fatjet FatJet_;
  bool are_GenJetAK8_loaded_;
  Genjetak8 GenJetAK8_;
  bool are_GenJet_loaded_;
  Genjet GenJet_;
  bool are_GenPart_loaded_;
  Genpart GenPart_;
  bool are_GenVisTau_loaded_;
  Genvistau GenVisTau_;
  bool are_LHEPdfWeight_loaded_;
  Lhepdfweight LHEPdfWeight_;
  bool are_LHEScaleWeight_loaded_;
  Lhescaleweight LHEScaleWeight_;
  bool are_Jet_loaded_;
  Jet Jet_;
  bool are_Muon_loaded_;
  Muon Muon_;
  bool are_Photon_loaded_;
  Photon Photon_;
  bool are_GenDressedLepton_loaded_;
  Gendressedlepton GenDressedLepton_;
  bool are_SoftActivityJet_loaded_;
  Softactivityjet SoftActivityJet_;
  bool are_SubJet_loaded_;
  Subjet SubJet_;
  bool are_Tau_loaded_;
  Tau Tau_;
  bool are_TrigObj_loaded_;
  Trigobj TrigObj_;
  bool are_OtherPV_loaded_;
  Otherpv OtherPV_;
  bool are_SV_loaded_;
  Sv SV_;
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

