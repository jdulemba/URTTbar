#ifndef URStreamer_h
#define URStreamer_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include <vector>
using namespace std;

class Je: public TLorentzVector{
friend class URStreamer;
public:
//  Je(const Float_ &i_area_,const Float_ &i_btagCMVA_,const Float_ &i_btagCSVV2_,const Float_ &i_btagDeepB_,const Float_ &i_btagDeepC_,const Float_ &i_chEmEF_,const Float_ &i_chHEF_,const Float_ &i_neEmEF_,const Float_ &i_neHEF_,const Float_ &i_qgl_,const Float_ &i_rawFactor_,const Float_ &i_bReg_,const Int_ &i_electronIdx1_,const Int_ &i_electronIdx2_,const Int_ &i_jetId_,const Int_ &i_muonIdx1_,const Int_ &i_muonIdx2_,const Int_ &i_nConstituents_,const Int_ &i_nElectrons_,const Int_ &i_nMuons_,const Int_ &i_puId_,const Int_ &i_genJetIdx_,const Int_ &i_hadronFlavour_,const Int_ &i_partonFlavour_,const UChar_ &i_cleanmask_):
//    
//  {}
  Je():
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
  Float_ area() const {return area_;}
  Float_ btagCMVA() const {return btagCMVA_;}
  Float_ btagCSVV2() const {return btagCSVV2_;}
  Float_ btagDeepB() const {return btagDeepB_;}
  Float_ btagDeepC() const {return btagDeepC_;}
  Float_ chEmEF() const {return chEmEF_;}
  Float_ chHEF() const {return chHEF_;}
  Float_ neEmEF() const {return neEmEF_;}
  Float_ neHEF() const {return neHEF_;}
  Float_ qgl() const {return qgl_;}
  Float_ rawFactor() const {return rawFactor_;}
  Float_ bReg() const {return bReg_;}
  Int_ electronIdx1() const {return electronIdx1_;}
  Int_ electronIdx2() const {return electronIdx2_;}
  Int_ jetId() const {return jetId_;}
  Int_ muonIdx1() const {return muonIdx1_;}
  Int_ muonIdx2() const {return muonIdx2_;}
  Int_ nConstituents() const {return nConstituents_;}
  Int_ nElectrons() const {return nElectrons_;}
  Int_ nMuons() const {return nMuons_;}
  Int_ puId() const {return puId_;}
  Int_ genJetIdx() const {return genJetIdx_;}
  Int_ hadronFlavour() const {return hadronFlavour_;}
  Int_ partonFlavour() const {return partonFlavour_;}
  UChar_ cleanmask() const {return cleanmask_;}
private:
  Float_ area_;
  Float_ btagCMVA_;
  Float_ btagCSVV2_;
  Float_ btagDeepB_;
  Float_ btagDeepC_;
  Float_ chEmEF_;
  Float_ chHEF_;
  Float_ neEmEF_;
  Float_ neHEF_;
  Float_ qgl_;
  Float_ rawFactor_;
  Float_ bReg_;
  Int_ electronIdx1_;
  Int_ electronIdx2_;
  Int_ jetId_;
  Int_ muonIdx1_;
  Int_ muonIdx2_;
  Int_ nConstituents_;
  Int_ nElectrons_;
  Int_ nMuons_;
  Int_ puId_;
  Int_ genJetIdx_;
  Int_ hadronFlavour_;
  Int_ partonFlavour_;
  UChar_ cleanmask_;
  void setarea(const Float_ value) {area_ = value;}
  void setbtagCMVA(const Float_ value) {btagCMVA_ = value;}
  void setbtagCSVV2(const Float_ value) {btagCSVV2_ = value;}
  void setbtagDeepB(const Float_ value) {btagDeepB_ = value;}
  void setbtagDeepC(const Float_ value) {btagDeepC_ = value;}
  void setchEmEF(const Float_ value) {chEmEF_ = value;}
  void setchHEF(const Float_ value) {chHEF_ = value;}
  void setneEmEF(const Float_ value) {neEmEF_ = value;}
  void setneHEF(const Float_ value) {neHEF_ = value;}
  void setqgl(const Float_ value) {qgl_ = value;}
  void setrawFactor(const Float_ value) {rawFactor_ = value;}
  void setbReg(const Float_ value) {bReg_ = value;}
  void setelectronIdx1(const Int_ value) {electronIdx1_ = value;}
  void setelectronIdx2(const Int_ value) {electronIdx2_ = value;}
  void setjetId(const Int_ value) {jetId_ = value;}
  void setmuonIdx1(const Int_ value) {muonIdx1_ = value;}
  void setmuonIdx2(const Int_ value) {muonIdx2_ = value;}
  void setnConstituents(const Int_ value) {nConstituents_ = value;}
  void setnElectrons(const Int_ value) {nElectrons_ = value;}
  void setnMuons(const Int_ value) {nMuons_ = value;}
  void setpuId(const Int_ value) {puId_ = value;}
  void setgenJetIdx(const Int_ value) {genJetIdx_ = value;}
  void sethadronFlavour(const Int_ value) {hadronFlavour_ = value;}
  void setpartonFlavour(const Int_ value) {partonFlavour_ = value;}
  void setcleanmask(const UChar_ value) {cleanmask_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Genjetak: public TLorentzVector{
friend class URStreamer;
public:
//  Genjetak(const Int_ &i_partonFlavour_,const UChar_ &i_hadronFlavour_):
//    
//  {}
  Genjetak():
    TLorentzVector(),
    partonFlavour_(0),
    hadronFlavour_(0)
  {}
  Int_ partonFlavour() const {return partonFlavour_;}
  UChar_ hadronFlavour() const {return hadronFlavour_;}
private:
  Int_ partonFlavour_;
  UChar_ hadronFlavour_;
  void setpartonFlavour(const Int_ value) {partonFlavour_ = value;}
  void sethadronFlavour(const UChar_ value) {hadronFlavour_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Genvista: public TLorentzVector{
friend class URStreamer;
public:
//  Genvista(const Int_ &i_charge_,const Int_ &i_genPartIdxMother_,const Int_ &i_status_):
//    
//  {}
  Genvista():
    TLorentzVector(),
    charge_(0),
    genPartIdxMother_(0),
    status_(0)
  {}
  Int_ charge() const {return charge_;}
  Int_ genPartIdxMother() const {return genPartIdxMother_;}
  Int_ status() const {return status_;}
private:
  Int_ charge_;
  Int_ genPartIdxMother_;
  Int_ status_;
  void setcharge(const Int_ value) {charge_ = value;}
  void setgenPartIdxMother(const Int_ value) {genPartIdxMother_ = value;}
  void setstatus(const Int_ value) {status_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Calomet: public TLorentzVector{
friend class URStreamer;
public:
//  Calomet(const Float_t &i_sumEt_):
//    
//  {}
  Calomet():
    TLorentzVector(),
    sumEt_(0)
  {}
  Float_t sumEt() const {return sumEt_;}
private:
  Float_t sumEt_;
  void setsumEt(const Float_t value) {sumEt_ = value;}
  void setLotentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Gendressedlepto: public TLorentzVector{
friend class URStreamer;
public:
//  Gendressedlepto(const Int_ &i_pdgId_):
//    
//  {}
  Gendressedlepto():
    TLorentzVector(),
    pdgId_(0)
  {}
  Int_ pdgId() const {return pdgId_;}
private:
  Int_ pdgId_;
  void setpdgId(const Int_ value) {pdgId_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Pv{
friend class URStreamer;
public:
//  Pv(const Float_t &i_ndof_,const Float_t &i_x_,const Float_t &i_y_,const Float_t &i_z_,const Float_t &i_chi2_,const Float_t &i_score_,const Int_t &i_npvs_,const Int_t &i_npvsGood_):
//    
//  {}
  Pv():
    ndof_(0),
    x_(0),
    y_(0),
    z_(0),
    chi2_(0),
    score_(0),
    npvs_(0),
    npvsGood_(0)
  {}
  Float_t ndof() const {return ndof_;}
  Float_t x() const {return x_;}
  Float_t y() const {return y_;}
  Float_t z() const {return z_;}
  Float_t chi2() const {return chi2_;}
  Float_t score() const {return score_;}
  Int_t npvs() const {return npvs_;}
  Int_t npvsGood() const {return npvsGood_;}
private:
  Float_t ndof_;
  Float_t x_;
  Float_t y_;
  Float_t z_;
  Float_t chi2_;
  Float_t score_;
  Int_t npvs_;
  Int_t npvsGood_;
  void setndof(const Float_t value) {ndof_ = value;}
  void setx(const Float_t value) {x_ = value;}
  void sety(const Float_t value) {y_ = value;}
  void setz(const Float_t value) {z_ = value;}
  void setchi2(const Float_t value) {chi2_ = value;}
  void setscore(const Float_t value) {score_ = value;}
  void setnpvs(const Int_t value) {npvs_ = value;}
  void setnpvsGood(const Int_t value) {npvsGood_ = value;}
};

class Generator{
friend class URStreamer;
public:
//  Generator(const Float_t &i_binvar_,const Float_t &i_scalePDF_,const Float_t &i_weight_,const Float_t &i_x1_,const Float_t &i_x2_,const Float_t &i_xpdf1_,const Float_t &i_xpdf2_,const Int_t &i_id1_,const Int_t &i_id2_):
//    
//  {}
  Generator():
    binvar_(0),
    scalePDF_(0),
    weight_(0),
    x1_(0),
    x2_(0),
    xpdf1_(0),
    xpdf2_(0),
    id1_(0),
    id2_(0)
  {}
  Float_t binvar() const {return binvar_;}
  Float_t scalePDF() const {return scalePDF_;}
  Float_t weight() const {return weight_;}
  Float_t x1() const {return x1_;}
  Float_t x2() const {return x2_;}
  Float_t xpdf1() const {return xpdf1_;}
  Float_t xpdf2() const {return xpdf2_;}
  Int_t id1() const {return id1_;}
  Int_t id2() const {return id2_;}
private:
  Float_t binvar_;
  Float_t scalePDF_;
  Float_t weight_;
  Float_t x1_;
  Float_t x2_;
  Float_t xpdf1_;
  Float_t xpdf2_;
  Int_t id1_;
  Int_t id2_;
  void setbinvar(const Float_t value) {binvar_ = value;}
  void setscalePDF(const Float_t value) {scalePDF_ = value;}
  void setweight(const Float_t value) {weight_ = value;}
  void setx1(const Float_t value) {x1_ = value;}
  void setx2(const Float_t value) {x2_ = value;}
  void setxpdf1(const Float_t value) {xpdf1_ = value;}
  void setxpdf2(const Float_t value) {xpdf2_ = value;}
  void setid1(const Int_t value) {id1_ = value;}
  void setid2(const Int_t value) {id2_ = value;}
};

class Trigob: public TLorentzVector{
friend class URStreamer;
public:
//  Trigob(const Float_ &i_l1pt_,const Float_ &i_l1pt_,const Float_ &i_l2pt_,const Int_ &i_id_,const Int_ &i_l1iso_,const Int_ &i_l1charge_,const Int_ &i_filterBits_):
//    
//  {}
  Trigob():
    TLorentzVector(),
    l1pt_(0),
    l1pt_(0),
    l2pt_(0),
    id_(0),
    l1iso_(0),
    l1charge_(0),
    filterBits_(0)
  {}
  Float_ l1pt() const {return l1pt_;}
  Float_ l1pt() const {return l1pt_;}
  Float_ l2pt() const {return l2pt_;}
  Int_ id() const {return id_;}
  Int_ l1iso() const {return l1iso_;}
  Int_ l1charge() const {return l1charge_;}
  Int_ filterBits() const {return filterBits_;}
private:
  Float_ l1pt_;
  Float_ l1pt_;
  Float_ l2pt_;
  Int_ id_;
  Int_ l1iso_;
  Int_ l1charge_;
  Int_ filterBits_;
  void setl1pt(const Float_ value) {l1pt_ = value;}
  void setl1pt(const Float_ value) {l1pt_ = value;}
  void setl2pt(const Float_ value) {l2pt_ = value;}
  void setid(const Int_ value) {id_ = value;}
  void setl1iso(const Int_ value) {l1iso_ = value;}
  void setl1charge(const Int_ value) {l1charge_ = value;}
  void setfilterBits(const Int_ value) {filterBits_ = value;}
  void setLotentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Photo: public TLorentzVector{
friend class URStreamer;
public:
//  Photo(const Float_ &i_eCorr_,const Float_ &i_energyErr_,const Float_ &i_hoe_,const Float_ &i_mvaID_,const Float_ &i_pfRelIso03_,const Float_ &i_pfRelIso03_,const Float_ &i_r9_,const Float_ &i_sieie_,const Int_ &i_charge_,const Int_ &i_cutBased_,const Int_ &i_electronIdx_,const Int_ &i_jetIdx_,const Int_ &i_pdgId_,const Int_ &i_vidNestedWPBitmap_,const Bool_ &i_electronVeto_,const Bool_ &i_mvaID_,const Bool_ &i_mvaID_,const Bool_ &i_pixelSeed_,const Int_ &i_genPartIdx_,const UChar_ &i_genPartFlav_,const UChar_ &i_cleanmask_):
//    
//  {}
  Photo():
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
  Float_ eCorr() const {return eCorr_;}
  Float_ energyErr() const {return energyErr_;}
  Float_ hoe() const {return hoe_;}
  Float_ mvaID() const {return mvaID_;}
  Float_ pfRelIso03() const {return pfRelIso03_;}
  Float_ pfRelIso03() const {return pfRelIso03_;}
  Float_ r9() const {return r9_;}
  Float_ sieie() const {return sieie_;}
  Int_ charge() const {return charge_;}
  Int_ cutBased() const {return cutBased_;}
  Int_ electronIdx() const {return electronIdx_;}
  Int_ jetIdx() const {return jetIdx_;}
  Int_ pdgId() const {return pdgId_;}
  Int_ vidNestedWPBitmap() const {return vidNestedWPBitmap_;}
  Bool_ electronVeto() const {return electronVeto_;}
  Bool_ mvaID() const {return mvaID_;}
  Bool_ mvaID() const {return mvaID_;}
  Bool_ pixelSeed() const {return pixelSeed_;}
  Int_ genPartIdx() const {return genPartIdx_;}
  UChar_ genPartFlav() const {return genPartFlav_;}
  UChar_ cleanmask() const {return cleanmask_;}
private:
  Float_ eCorr_;
  Float_ energyErr_;
  Float_ hoe_;
  Float_ mvaID_;
  Float_ pfRelIso03_;
  Float_ pfRelIso03_;
  Float_ r9_;
  Float_ sieie_;
  Int_ charge_;
  Int_ cutBased_;
  Int_ electronIdx_;
  Int_ jetIdx_;
  Int_ pdgId_;
  Int_ vidNestedWPBitmap_;
  Bool_ electronVeto_;
  Bool_ mvaID_;
  Bool_ mvaID_;
  Bool_ pixelSeed_;
  Int_ genPartIdx_;
  UChar_ genPartFlav_;
  UChar_ cleanmask_;
  void seteCorr(const Float_ value) {eCorr_ = value;}
  void setenergyErr(const Float_ value) {energyErr_ = value;}
  void sethoe(const Float_ value) {hoe_ = value;}
  void setmvaID(const Float_ value) {mvaID_ = value;}
  void setpfRelIso03(const Float_ value) {pfRelIso03_ = value;}
  void setpfRelIso03(const Float_ value) {pfRelIso03_ = value;}
  void setr9(const Float_ value) {r9_ = value;}
  void setsieie(const Float_ value) {sieie_ = value;}
  void setcharge(const Int_ value) {charge_ = value;}
  void setcutBased(const Int_ value) {cutBased_ = value;}
  void setelectronIdx(const Int_ value) {electronIdx_ = value;}
  void setjetIdx(const Int_ value) {jetIdx_ = value;}
  void setpdgId(const Int_ value) {pdgId_ = value;}
  void setvidNestedWPBitmap(const Int_ value) {vidNestedWPBitmap_ = value;}
  void setelectronVeto(const Bool_ value) {electronVeto_ = value;}
  void setmvaID(const Bool_ value) {mvaID_ = value;}
  void setmvaID(const Bool_ value) {mvaID_ = value;}
  void setpixelSeed(const Bool_ value) {pixelSeed_ = value;}
  void setgenPartIdx(const Int_ value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_ value) {genPartFlav_ = value;}
  void setcleanmask(const UChar_ value) {cleanmask_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Genje: public TLorentzVector{
friend class URStreamer;
public:
//  Genje(const Int_ &i_partonFlavour_,const UChar_ &i_hadronFlavour_,const Int_ &i_partonFlavour_,const UChar_ &i_hadronFlavour_):
//    
//  {}
  Genje():
    TLorentzVector(),
    partonFlavour_(0),
    hadronFlavour_(0),
    partonFlavour_(0),
    hadronFlavour_(0)
  {}
  Int_ partonFlavour() const {return partonFlavour_;}
  UChar_ hadronFlavour() const {return hadronFlavour_;}
  Int_ partonFlavour() const {return partonFlavour_;}
  UChar_ hadronFlavour() const {return hadronFlavour_;}
private:
  Int_ partonFlavour_;
  UChar_ hadronFlavour_;
  Int_ partonFlavour_;
  UChar_ hadronFlavour_;
  void setpartonFlavour(const Int_ value) {partonFlavour_ = value;}
  void sethadronFlavour(const UChar_ value) {hadronFlavour_ = value;}
  void setpartonFlavour(const Int_ value) {partonFlavour_ = value;}
  void sethadronFlavour(const UChar_ value) {hadronFlavour_ = value;}
};

class Rawmet: public TLorentzVector{
friend class URStreamer;
public:
//  Rawmet(const Float_t &i_sumEt_):
//    
//  {}
  Rawmet():
    TLorentzVector(),
    sumEt_(0)
  {}
  Float_t sumEt() const {return sumEt_;}
private:
  Float_t sumEt_;
  void setsumEt(const Float_t value) {sumEt_ = value;}
  void setLotentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Electro: public TLorentzVector{
friend class URStreamer;
public:
//  Electro(const Float_ &i_deltaEtaSC_,const Float_ &i_dr03EcalRecHitSumEt_,const Float_ &i_dr03HcalDepth1TowerSumEt_,const Float_ &i_dr03TkSumPt_,const Float_ &i_dxy_,const Float_ &i_dxyErr_,const Float_ &i_dz_,const Float_ &i_dzErr_,const Float_ &i_eCorr_,const Float_ &i_eInvMinusPInv_,const Float_ &i_energyErr_,const Float_ &i_hoe_,const Float_ &i_ip3d_,const Float_ &i_miniPFRelIso_,const Float_ &i_miniPFRelIso_,const Float_ &i_mvaSpring16GP_,const Float_ &i_mvaSpring16HZZ_,const Float_ &i_pfRelIso03_,const Float_ &i_pfRelIso03_,const Float_ &i_r9_,const Float_ &i_sieie_,const Float_ &i_sip3d_,const Float_ &i_mvaTTH_,const Int_ &i_charge_,const Int_ &i_cutBased_,const Int_ &i_cutBased_,const Int_ &i_jetIdx_,const Int_ &i_pdgId_,const Int_ &i_photonIdx_,const Int_ &i_tightCharge_,const Int_ &i_vidNestedWPBitmap_,const Bool_ &i_convVeto_,const Bool_ &i_cutBased_,const Bool_ &i_isPFcand_,const UChar_ &i_lostHits_,const Bool_ &i_mvaSpring16GP_,const Bool_ &i_mvaSpring16GP_,const Bool_ &i_mvaSpring16HZZ_,const Int_ &i_genPartIdx_,const UChar_ &i_genPartFlav_,const UChar_ &i_cleanmask_):
//    
//  {}
  Electro():
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
  Float_ deltaEtaSC() const {return deltaEtaSC_;}
  Float_ dr03EcalRecHitSumEt() const {return dr03EcalRecHitSumEt_;}
  Float_ dr03HcalDepth1TowerSumEt() const {return dr03HcalDepth1TowerSumEt_;}
  Float_ dr03TkSumPt() const {return dr03TkSumPt_;}
  Float_ dxy() const {return dxy_;}
  Float_ dxyErr() const {return dxyErr_;}
  Float_ dz() const {return dz_;}
  Float_ dzErr() const {return dzErr_;}
  Float_ eCorr() const {return eCorr_;}
  Float_ eInvMinusPInv() const {return eInvMinusPInv_;}
  Float_ energyErr() const {return energyErr_;}
  Float_ hoe() const {return hoe_;}
  Float_ ip3d() const {return ip3d_;}
  Float_ miniPFRelIso() const {return miniPFRelIso_;}
  Float_ miniPFRelIso() const {return miniPFRelIso_;}
  Float_ mvaSpring16GP() const {return mvaSpring16GP_;}
  Float_ mvaSpring16HZZ() const {return mvaSpring16HZZ_;}
  Float_ pfRelIso03() const {return pfRelIso03_;}
  Float_ pfRelIso03() const {return pfRelIso03_;}
  Float_ r9() const {return r9_;}
  Float_ sieie() const {return sieie_;}
  Float_ sip3d() const {return sip3d_;}
  Float_ mvaTTH() const {return mvaTTH_;}
  Int_ charge() const {return charge_;}
  Int_ cutBased() const {return cutBased_;}
  Int_ cutBased() const {return cutBased_;}
  Int_ jetIdx() const {return jetIdx_;}
  Int_ pdgId() const {return pdgId_;}
  Int_ photonIdx() const {return photonIdx_;}
  Int_ tightCharge() const {return tightCharge_;}
  Int_ vidNestedWPBitmap() const {return vidNestedWPBitmap_;}
  Bool_ convVeto() const {return convVeto_;}
  Bool_ cutBased() const {return cutBased_;}
  Bool_ isPFcand() const {return isPFcand_;}
  UChar_ lostHits() const {return lostHits_;}
  Bool_ mvaSpring16GP() const {return mvaSpring16GP_;}
  Bool_ mvaSpring16GP() const {return mvaSpring16GP_;}
  Bool_ mvaSpring16HZZ() const {return mvaSpring16HZZ_;}
  Int_ genPartIdx() const {return genPartIdx_;}
  UChar_ genPartFlav() const {return genPartFlav_;}
  UChar_ cleanmask() const {return cleanmask_;}
private:
  Float_ deltaEtaSC_;
  Float_ dr03EcalRecHitSumEt_;
  Float_ dr03HcalDepth1TowerSumEt_;
  Float_ dr03TkSumPt_;
  Float_ dxy_;
  Float_ dxyErr_;
  Float_ dz_;
  Float_ dzErr_;
  Float_ eCorr_;
  Float_ eInvMinusPInv_;
  Float_ energyErr_;
  Float_ hoe_;
  Float_ ip3d_;
  Float_ miniPFRelIso_;
  Float_ miniPFRelIso_;
  Float_ mvaSpring16GP_;
  Float_ mvaSpring16HZZ_;
  Float_ pfRelIso03_;
  Float_ pfRelIso03_;
  Float_ r9_;
  Float_ sieie_;
  Float_ sip3d_;
  Float_ mvaTTH_;
  Int_ charge_;
  Int_ cutBased_;
  Int_ cutBased_;
  Int_ jetIdx_;
  Int_ pdgId_;
  Int_ photonIdx_;
  Int_ tightCharge_;
  Int_ vidNestedWPBitmap_;
  Bool_ convVeto_;
  Bool_ cutBased_;
  Bool_ isPFcand_;
  UChar_ lostHits_;
  Bool_ mvaSpring16GP_;
  Bool_ mvaSpring16GP_;
  Bool_ mvaSpring16HZZ_;
  Int_ genPartIdx_;
  UChar_ genPartFlav_;
  UChar_ cleanmask_;
  void setdeltaEtaSC(const Float_ value) {deltaEtaSC_ = value;}
  void setdr03EcalRecHitSumEt(const Float_ value) {dr03EcalRecHitSumEt_ = value;}
  void setdr03HcalDepth1TowerSumEt(const Float_ value) {dr03HcalDepth1TowerSumEt_ = value;}
  void setdr03TkSumPt(const Float_ value) {dr03TkSumPt_ = value;}
  void setdxy(const Float_ value) {dxy_ = value;}
  void setdxyErr(const Float_ value) {dxyErr_ = value;}
  void setdz(const Float_ value) {dz_ = value;}
  void setdzErr(const Float_ value) {dzErr_ = value;}
  void seteCorr(const Float_ value) {eCorr_ = value;}
  void seteInvMinusPInv(const Float_ value) {eInvMinusPInv_ = value;}
  void setenergyErr(const Float_ value) {energyErr_ = value;}
  void sethoe(const Float_ value) {hoe_ = value;}
  void setip3d(const Float_ value) {ip3d_ = value;}
  void setminiPFRelIso(const Float_ value) {miniPFRelIso_ = value;}
  void setminiPFRelIso(const Float_ value) {miniPFRelIso_ = value;}
  void setmvaSpring16GP(const Float_ value) {mvaSpring16GP_ = value;}
  void setmvaSpring16HZZ(const Float_ value) {mvaSpring16HZZ_ = value;}
  void setpfRelIso03(const Float_ value) {pfRelIso03_ = value;}
  void setpfRelIso03(const Float_ value) {pfRelIso03_ = value;}
  void setr9(const Float_ value) {r9_ = value;}
  void setsieie(const Float_ value) {sieie_ = value;}
  void setsip3d(const Float_ value) {sip3d_ = value;}
  void setmvaTTH(const Float_ value) {mvaTTH_ = value;}
  void setcharge(const Int_ value) {charge_ = value;}
  void setcutBased(const Int_ value) {cutBased_ = value;}
  void setcutBased(const Int_ value) {cutBased_ = value;}
  void setjetIdx(const Int_ value) {jetIdx_ = value;}
  void setpdgId(const Int_ value) {pdgId_ = value;}
  void setphotonIdx(const Int_ value) {photonIdx_ = value;}
  void settightCharge(const Int_ value) {tightCharge_ = value;}
  void setvidNestedWPBitmap(const Int_ value) {vidNestedWPBitmap_ = value;}
  void setconvVeto(const Bool_ value) {convVeto_ = value;}
  void setcutBased(const Bool_ value) {cutBased_ = value;}
  void setisPFcand(const Bool_ value) {isPFcand_ = value;}
  void setlostHits(const UChar_ value) {lostHits_ = value;}
  void setmvaSpring16GP(const Bool_ value) {mvaSpring16GP_ = value;}
  void setmvaSpring16GP(const Bool_ value) {mvaSpring16GP_ = value;}
  void setmvaSpring16HZZ(const Bool_ value) {mvaSpring16HZZ_ = value;}
  void setgenPartIdx(const Int_ value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_ value) {genPartFlav_ = value;}
  void setcleanmask(const UChar_ value) {cleanmask_ = value;}
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

class L1Simulation{
friend class URStreamer;
public:
//  L1Simulation(const Bool_t &i_step_):
//    
//  {}
  L1Simulation():
    step_(0)
  {}
  Bool_t step() const {return step_;}
private:
  Bool_t step_;
  void setstep(const Bool_t value) {step_ = value;}
};

class Genpar: public TLorentzVector{
friend class URStreamer;
public:
//  Genpar(const Int_ &i_genPartIdxMother_,const Int_ &i_pdgId_,const Int_ &i_status_,const Int_ &i_statusFlags_):
//    
//  {}
  Genpar():
    TLorentzVector(),
    genPartIdxMother_(0),
    pdgId_(0),
    status_(0),
    statusFlags_(0)
  {}
  Int_ genPartIdxMother() const {return genPartIdxMother_;}
  Int_ pdgId() const {return pdgId_;}
  Int_ status() const {return status_;}
  Int_ statusFlags() const {return statusFlags_;}
private:
  Int_ genPartIdxMother_;
  Int_ pdgId_;
  Int_ status_;
  Int_ statusFlags_;
  void setgenPartIdxMother(const Int_ value) {genPartIdxMother_ = value;}
  void setpdgId(const Int_ value) {pdgId_ = value;}
  void setstatus(const Int_ value) {status_ = value;}
  void setstatusFlags(const Int_ value) {statusFlags_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Lhe{
friend class URStreamer;
public:
//  Lhe(const Float_t &i_originalXWGTUP_,const Float_ &i_LHEPdfWeight_,const Float_ &i_LHEScaleWeight_,const Float_t &i_HT_,const Float_t &i_HTIncoming_,const Float_t &i_Vpt_,const UChar_t &i_Njets_,const UChar_t &i_Nb_,const UChar_t &i_Nc_,const UChar_t &i_Nuds_,const UChar_t &i_Nglu_,const UChar_t &i_NpNLO_,const UChar_t &i_NpLO_):
//    
//  {}
  Lhe():
    originalXWGTUP_(0),
    LHEPdfWeight_(0),
    LHEScaleWeight_(0),
    HT_(0),
    HTIncoming_(0),
    Vpt_(0),
    Njets_(0),
    Nb_(0),
    Nc_(0),
    Nuds_(0),
    Nglu_(0),
    NpNLO_(0),
    NpLO_(0)
  {}
  Float_t originalXWGTUP() const {return originalXWGTUP_;}
  Float_ LHEPdfWeight() const {return LHEPdfWeight_;}
  Float_ LHEScaleWeight() const {return LHEScaleWeight_;}
  Float_t HT() const {return HT_;}
  Float_t HTIncoming() const {return HTIncoming_;}
  Float_t Vpt() const {return Vpt_;}
  UChar_t Njets() const {return Njets_;}
  UChar_t Nb() const {return Nb_;}
  UChar_t Nc() const {return Nc_;}
  UChar_t Nuds() const {return Nuds_;}
  UChar_t Nglu() const {return Nglu_;}
  UChar_t NpNLO() const {return NpNLO_;}
  UChar_t NpLO() const {return NpLO_;}
private:
  Float_t originalXWGTUP_;
  Float_ LHEPdfWeight_;
  Float_ LHEScaleWeight_;
  Float_t HT_;
  Float_t HTIncoming_;
  Float_t Vpt_;
  UChar_t Njets_;
  UChar_t Nb_;
  UChar_t Nc_;
  UChar_t Nuds_;
  UChar_t Nglu_;
  UChar_t NpNLO_;
  UChar_t NpLO_;
  void setoriginalXWGTUP(const Float_t value) {originalXWGTUP_ = value;}
  void setLHEPdfWeight(const Float_ value) {LHEPdfWeight_ = value;}
  void setLHEScaleWeight(const Float_ value) {LHEScaleWeight_ = value;}
  void setHT(const Float_t value) {HT_ = value;}
  void setHTIncoming(const Float_t value) {HTIncoming_ = value;}
  void setVpt(const Float_t value) {Vpt_ = value;}
  void setNjets(const UChar_t value) {Njets_ = value;}
  void setNb(const UChar_t value) {Nb_ = value;}
  void setNc(const UChar_t value) {Nc_ = value;}
  void setNuds(const UChar_t value) {Nuds_ = value;}
  void setNglu(const UChar_t value) {Nglu_ = value;}
  void setNpNLO(const UChar_t value) {NpNLO_ = value;}
  void setNpLO(const UChar_t value) {NpLO_ = value;}
};

class Tkmet: public TLorentzVector{
friend class URStreamer;
public:
//  Tkmet(const Float_t &i_sumEt_):
//    
//  {}
  Tkmet():
    TLorentzVector(),
    sumEt_(0)
  {}
  Float_t sumEt() const {return sumEt_;}
private:
  Float_t sumEt_;
  void setsumEt(const Float_t value) {sumEt_ = value;}
  void setLotentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Ta: public TLorentzVector{
friend class URStreamer;
public:
//  Ta(const Float_ &i_chargedIso_,const Float_ &i_dxy_,const Float_ &i_dz_,const Float_ &i_footprintCorr_,const Float_ &i_leadTkDeltaEta_,const Float_ &i_leadTkDeltaPhi_,const Float_ &i_leadTkPtOverTauPt_,const Float_ &i_neutralIso_,const Float_ &i_photonsOutsideSignalCone_,const Float_ &i_puCorr_,const Float_ &i_rawAntiEle_,const Float_ &i_rawIso_,const Float_ &i_rawMVAnewDM_,const Float_ &i_rawMVAoldDM_,const Float_ &i_rawMVAoldDMdR03_,const Int_ &i_charge_,const Int_ &i_decayMode_,const Int_ &i_jetIdx_,const Int_ &i_rawAntiEleCat_,const UChar_ &i_idAntiEle_,const UChar_ &i_idAntiMu_,const Bool_ &i_idDecayMode_,const Bool_ &i_idDecayModeNewDMs_,const UChar_ &i_idMVAnewDM_,const UChar_ &i_idMVAoldDM_,const UChar_ &i_idMVAoldDMdR03_,const UChar_ &i_cleanmask_,const Int_ &i_genPartIdx_,const UChar_ &i_genPartFlav_):
//    
//  {}
  Ta():
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
  Float_ chargedIso() const {return chargedIso_;}
  Float_ dxy() const {return dxy_;}
  Float_ dz() const {return dz_;}
  Float_ footprintCorr() const {return footprintCorr_;}
  Float_ leadTkDeltaEta() const {return leadTkDeltaEta_;}
  Float_ leadTkDeltaPhi() const {return leadTkDeltaPhi_;}
  Float_ leadTkPtOverTauPt() const {return leadTkPtOverTauPt_;}
  Float_ neutralIso() const {return neutralIso_;}
  Float_ photonsOutsideSignalCone() const {return photonsOutsideSignalCone_;}
  Float_ puCorr() const {return puCorr_;}
  Float_ rawAntiEle() const {return rawAntiEle_;}
  Float_ rawIso() const {return rawIso_;}
  Float_ rawMVAnewDM() const {return rawMVAnewDM_;}
  Float_ rawMVAoldDM() const {return rawMVAoldDM_;}
  Float_ rawMVAoldDMdR03() const {return rawMVAoldDMdR03_;}
  Int_ charge() const {return charge_;}
  Int_ decayMode() const {return decayMode_;}
  Int_ jetIdx() const {return jetIdx_;}
  Int_ rawAntiEleCat() const {return rawAntiEleCat_;}
  UChar_ idAntiEle() const {return idAntiEle_;}
  UChar_ idAntiMu() const {return idAntiMu_;}
  Bool_ idDecayMode() const {return idDecayMode_;}
  Bool_ idDecayModeNewDMs() const {return idDecayModeNewDMs_;}
  UChar_ idMVAnewDM() const {return idMVAnewDM_;}
  UChar_ idMVAoldDM() const {return idMVAoldDM_;}
  UChar_ idMVAoldDMdR03() const {return idMVAoldDMdR03_;}
  UChar_ cleanmask() const {return cleanmask_;}
  Int_ genPartIdx() const {return genPartIdx_;}
  UChar_ genPartFlav() const {return genPartFlav_;}
private:
  Float_ chargedIso_;
  Float_ dxy_;
  Float_ dz_;
  Float_ footprintCorr_;
  Float_ leadTkDeltaEta_;
  Float_ leadTkDeltaPhi_;
  Float_ leadTkPtOverTauPt_;
  Float_ neutralIso_;
  Float_ photonsOutsideSignalCone_;
  Float_ puCorr_;
  Float_ rawAntiEle_;
  Float_ rawIso_;
  Float_ rawMVAnewDM_;
  Float_ rawMVAoldDM_;
  Float_ rawMVAoldDMdR03_;
  Int_ charge_;
  Int_ decayMode_;
  Int_ jetIdx_;
  Int_ rawAntiEleCat_;
  UChar_ idAntiEle_;
  UChar_ idAntiMu_;
  Bool_ idDecayMode_;
  Bool_ idDecayModeNewDMs_;
  UChar_ idMVAnewDM_;
  UChar_ idMVAoldDM_;
  UChar_ idMVAoldDMdR03_;
  UChar_ cleanmask_;
  Int_ genPartIdx_;
  UChar_ genPartFlav_;
  void setchargedIso(const Float_ value) {chargedIso_ = value;}
  void setdxy(const Float_ value) {dxy_ = value;}
  void setdz(const Float_ value) {dz_ = value;}
  void setfootprintCorr(const Float_ value) {footprintCorr_ = value;}
  void setleadTkDeltaEta(const Float_ value) {leadTkDeltaEta_ = value;}
  void setleadTkDeltaPhi(const Float_ value) {leadTkDeltaPhi_ = value;}
  void setleadTkPtOverTauPt(const Float_ value) {leadTkPtOverTauPt_ = value;}
  void setneutralIso(const Float_ value) {neutralIso_ = value;}
  void setphotonsOutsideSignalCone(const Float_ value) {photonsOutsideSignalCone_ = value;}
  void setpuCorr(const Float_ value) {puCorr_ = value;}
  void setrawAntiEle(const Float_ value) {rawAntiEle_ = value;}
  void setrawIso(const Float_ value) {rawIso_ = value;}
  void setrawMVAnewDM(const Float_ value) {rawMVAnewDM_ = value;}
  void setrawMVAoldDM(const Float_ value) {rawMVAoldDM_ = value;}
  void setrawMVAoldDMdR03(const Float_ value) {rawMVAoldDMdR03_ = value;}
  void setcharge(const Int_ value) {charge_ = value;}
  void setdecayMode(const Int_ value) {decayMode_ = value;}
  void setjetIdx(const Int_ value) {jetIdx_ = value;}
  void setrawAntiEleCat(const Int_ value) {rawAntiEleCat_ = value;}
  void setidAntiEle(const UChar_ value) {idAntiEle_ = value;}
  void setidAntiMu(const UChar_ value) {idAntiMu_ = value;}
  void setidDecayMode(const Bool_ value) {idDecayMode_ = value;}
  void setidDecayModeNewDMs(const Bool_ value) {idDecayModeNewDMs_ = value;}
  void setidMVAnewDM(const UChar_ value) {idMVAnewDM_ = value;}
  void setidMVAoldDM(const UChar_ value) {idMVAoldDM_ = value;}
  void setidMVAoldDMdR03(const UChar_ value) {idMVAoldDMdR03_ = value;}
  void setcleanmask(const UChar_ value) {cleanmask_ = value;}
  void setgenPartIdx(const Int_ value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_ value) {genPartFlav_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Puppimet: public TLorentzVector{
friend class URStreamer;
public:
//  Puppimet(const Float_t &i_sumEt_):
//    
//  {}
  Puppimet():
    TLorentzVector(),
    sumEt_(0)
  {}
  Float_t sumEt() const {return sumEt_;}
private:
  Float_t sumEt_;
  void setsumEt(const Float_t value) {sumEt_ = value;}
  void setLotentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Muo: public TLorentzVector{
friend class URStreamer;
public:
//  Muo(const Float_ &i_dxy_,const Float_ &i_dxyErr_,const Float_ &i_dz_,const Float_ &i_dzErr_,const Float_ &i_ip3d_,const Float_ &i_miniPFRelIso_,const Float_ &i_miniPFRelIso_,const Float_ &i_pfRelIso03_,const Float_ &i_pfRelIso03_,const Float_ &i_pfRelIso04_,const Float_ &i_ptErr_,const Float_ &i_segmentComp_,const Float_ &i_sip3d_,const Float_ &i_mvaTTH_,const Int_ &i_charge_,const Int_ &i_jetIdx_,const Int_ &i_nStations_,const Int_ &i_nTrackerLayers_,const Int_ &i_pdgId_,const Int_ &i_tightCharge_,const UChar_ &i_highPtId_,const Bool_ &i_isPFcand_,const Bool_ &i_mediumId_,const Bool_ &i_softId_,const Bool_ &i_tightId_,const Int_ &i_genPartIdx_,const UChar_ &i_genPartFlav_,const UChar_ &i_cleanmask_):
//    
//  {}
  Muo():
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
  Float_ dxy() const {return dxy_;}
  Float_ dxyErr() const {return dxyErr_;}
  Float_ dz() const {return dz_;}
  Float_ dzErr() const {return dzErr_;}
  Float_ ip3d() const {return ip3d_;}
  Float_ miniPFRelIso() const {return miniPFRelIso_;}
  Float_ miniPFRelIso() const {return miniPFRelIso_;}
  Float_ pfRelIso03() const {return pfRelIso03_;}
  Float_ pfRelIso03() const {return pfRelIso03_;}
  Float_ pfRelIso04() const {return pfRelIso04_;}
  Float_ ptErr() const {return ptErr_;}
  Float_ segmentComp() const {return segmentComp_;}
  Float_ sip3d() const {return sip3d_;}
  Float_ mvaTTH() const {return mvaTTH_;}
  Int_ charge() const {return charge_;}
  Int_ jetIdx() const {return jetIdx_;}
  Int_ nStations() const {return nStations_;}
  Int_ nTrackerLayers() const {return nTrackerLayers_;}
  Int_ pdgId() const {return pdgId_;}
  Int_ tightCharge() const {return tightCharge_;}
  UChar_ highPtId() const {return highPtId_;}
  Bool_ isPFcand() const {return isPFcand_;}
  Bool_ mediumId() const {return mediumId_;}
  Bool_ softId() const {return softId_;}
  Bool_ tightId() const {return tightId_;}
  Int_ genPartIdx() const {return genPartIdx_;}
  UChar_ genPartFlav() const {return genPartFlav_;}
  UChar_ cleanmask() const {return cleanmask_;}
private:
  Float_ dxy_;
  Float_ dxyErr_;
  Float_ dz_;
  Float_ dzErr_;
  Float_ ip3d_;
  Float_ miniPFRelIso_;
  Float_ miniPFRelIso_;
  Float_ pfRelIso03_;
  Float_ pfRelIso03_;
  Float_ pfRelIso04_;
  Float_ ptErr_;
  Float_ segmentComp_;
  Float_ sip3d_;
  Float_ mvaTTH_;
  Int_ charge_;
  Int_ jetIdx_;
  Int_ nStations_;
  Int_ nTrackerLayers_;
  Int_ pdgId_;
  Int_ tightCharge_;
  UChar_ highPtId_;
  Bool_ isPFcand_;
  Bool_ mediumId_;
  Bool_ softId_;
  Bool_ tightId_;
  Int_ genPartIdx_;
  UChar_ genPartFlav_;
  UChar_ cleanmask_;
  void setdxy(const Float_ value) {dxy_ = value;}
  void setdxyErr(const Float_ value) {dxyErr_ = value;}
  void setdz(const Float_ value) {dz_ = value;}
  void setdzErr(const Float_ value) {dzErr_ = value;}
  void setip3d(const Float_ value) {ip3d_ = value;}
  void setminiPFRelIso(const Float_ value) {miniPFRelIso_ = value;}
  void setminiPFRelIso(const Float_ value) {miniPFRelIso_ = value;}
  void setpfRelIso03(const Float_ value) {pfRelIso03_ = value;}
  void setpfRelIso03(const Float_ value) {pfRelIso03_ = value;}
  void setpfRelIso04(const Float_ value) {pfRelIso04_ = value;}
  void setptErr(const Float_ value) {ptErr_ = value;}
  void setsegmentComp(const Float_ value) {segmentComp_ = value;}
  void setsip3d(const Float_ value) {sip3d_ = value;}
  void setmvaTTH(const Float_ value) {mvaTTH_ = value;}
  void setcharge(const Int_ value) {charge_ = value;}
  void setjetIdx(const Int_ value) {jetIdx_ = value;}
  void setnStations(const Int_ value) {nStations_ = value;}
  void setnTrackerLayers(const Int_ value) {nTrackerLayers_ = value;}
  void setpdgId(const Int_ value) {pdgId_ = value;}
  void settightCharge(const Int_ value) {tightCharge_ = value;}
  void sethighPtId(const UChar_ value) {highPtId_ = value;}
  void setisPFcand(const Bool_ value) {isPFcand_ = value;}
  void setmediumId(const Bool_ value) {mediumId_ = value;}
  void setsoftId(const Bool_ value) {softId_ = value;}
  void settightId(const Bool_ value) {tightId_ = value;}
  void setgenPartIdx(const Int_ value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_ value) {genPartFlav_ = value;}
  void setcleanmask(const UChar_ value) {cleanmask_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Otherp{
friend class URStreamer;
public:
//  Otherp(const Float_ &i_z_):
//    
//  {}
  Otherp():
    z_(0)
  {}
  Float_ z() const {return z_;}
private:
  Float_ z_;
  void setz(const Float_ value) {z_ = value;}
};

class Hlt{
friend class URStreamer;
public:
//  Hlt(const Bool_t &i_HLTriggerFirstPath_,const Bool_t &i_AK8PFJet360_,const Bool_t &i_AK8PFJet400_,const Bool_t &i_AK8PFHT750_,const Bool_t &i_AK8PFHT800_,const Bool_t &i_AK8DiPFJet300_,const Bool_t &i_AK8DiPFJet280_,const Bool_t &i_AK8DiPFJet300_,const Bool_t &i_AK8DiPFJet300_,const Bool_t &i_AK8PFHT700_,const Bool_t &i_AK8PFHT650_,const Bool_t &i_AK8PFHT600_,const Bool_t &i_AK8DiPFJet280_,const Bool_t &i_AK8DiPFJet250_,const Bool_t &i_AK8DiPFJet280_,const Bool_t &i_AK8DiPFJet250_,const Bool_t &i_CaloJet260_,const Bool_t &i_CaloJet500_,const Bool_t &i_Dimuon13_,const Bool_t &i_Dimuon13_,const Bool_t &i_Dimuon20_,const Bool_t &i_DoubleEle24_,const Bool_t &i_DoubleEle25_,const Bool_t &i_DoubleEle33_,const Bool_t &i_DoubleEle33_,const Bool_t &i_DoubleEle33_,const Bool_t &i_DoubleEle33_,const Bool_t &i_DoubleMediumCombinedIsoPFTau35_,const Bool_t &i_DoubleTightCombinedIsoPFTau35_,const Bool_t &i_DoubleMediumCombinedIsoPFTau40_,const Bool_t &i_DoubleTightCombinedIsoPFTau40_,const Bool_t &i_DoubleMediumCombinedIsoPFTau40_,const Bool_t &i_DoubleTightCombinedIsoPFTau40_,const Bool_t &i_DoubleMediumIsoPFTau35_,const Bool_t &i_DoubleMediumIsoPFTau40_,const Bool_t &i_DoubleMediumIsoPFTau40_,const Bool_t &i_DoubleEle37_,const Bool_t &i_DoubleMu33NoFiltersNoVtx_,const Bool_t &i_DoubleMu38NoFiltersNoVtx_,const Bool_t &i_DoubleMu23NoFiltersNoVtxDisplaced_,const Bool_t &i_DoubleMu28NoFiltersNoVtxDisplaced_,const Bool_t &i_DoubleMu0_,const Bool_t &i_DoubleMu4_,const Bool_t &i_DoubleMu4_,const Bool_t &i_DoubleMu4_,const Bool_t &i_DoubleMu4_,const Bool_t &i_DoubleMu3_,const Bool_t &i_DoubleMu4_,const Bool_t &i_Mu7p5_,const Bool_t &i_Mu7p5_,const Bool_t &i_Mu7p5_,const Bool_t &i_Mu7p5_,const Bool_t &i_Mu7p5_,const Bool_t &i_Mu7p5_,const Bool_t &i_Mu7p5_,const Bool_t &i_Mu7p5_,const Bool_t &i_Dimuon0er16_,const Bool_t &i_Dimuon0er16_,const Bool_t &i_Dimuon6_,const Bool_t &i_Photon150_,const Bool_t &i_Photon90_,const Bool_t &i_HT250_,const Bool_t &i_DoublePhoton60_,const Bool_t &i_DoublePhoton85_,const Bool_t &i_Ele17_,const Bool_t &i_Ele20_,const Bool_t &i_Ele22_,const Bool_t &i_Ele22_,const Bool_t &i_Ele22_,const Bool_t &i_Ele23_,const Bool_t &i_Ele23_,const Bool_t &i_Ele24_,const Bool_t &i_Ele24_,const Bool_t &i_Ele24_,const Bool_t &i_Ele24_,const Bool_t &i_Ele25_,const Bool_t &i_Ele25_,const Bool_t &i_Ele25_,const Bool_t &i_Ele27_,const Bool_t &i_Ele27_,const Bool_t &i_Ele27_,const Bool_t &i_Ele27_,const Bool_t &i_Ele27_,const Bool_t &i_Ele27_,const Bool_t &i_Ele27_,const Bool_t &i_Ele30_,const Bool_t &i_Ele30_,const Bool_t &i_Ele30_,const Bool_t &i_Ele32_,const Bool_t &i_Ele32_,const Bool_t &i_Ele32_,const Bool_t &i_Ele32_,const Bool_t &i_Ele35_,const Bool_t &i_Ele35_,const Bool_t &i_Ele36_,const Bool_t &i_Ele45_,const Bool_t &i_Ele45_,const Bool_t &i_Ele45_,const Bool_t &i_Ele105_,const Bool_t &i_Ele30WP60_,const Bool_t &i_Ele30WP60_,const Bool_t &i_HT200_,const Bool_t &i_HT275_,const Bool_t &i_HT325_,const Bool_t &i_HT425_,const Bool_t &i_HT575_,const Bool_t &i_HT410to430_,const Bool_t &i_HT430to450_,const Bool_t &i_HT450to470_,const Bool_t &i_HT470to500_,const Bool_t &i_HT500to550_,const Bool_t &i_HT550to650_,const Bool_t &i_HT650_,const Bool_t &i_Mu16_,const Bool_t &i_IsoMu16_,const Bool_t &i_IsoMu16_,const Bool_t &i_IsoMu17_,const Bool_t &i_IsoMu17_,const Bool_t &i_IsoMu17_,const Bool_t &i_DoubleIsoMu17_,const Bool_t &i_DoubleIsoMu17_,const Bool_t &i_IsoMu18_,const Bool_t &i_IsoMu19_,const Bool_t &i_IsoMu19_,const Bool_t &i_IsoMu19_,const Bool_t &i_IsoMu19_,const Bool_t &i_IsoMu19_,const Bool_t &i_IsoMu19_,const Bool_t &i_IsoMu21_,const Bool_t &i_IsoMu21_,const Bool_t &i_IsoMu20_,const Bool_t &i_IsoMu21_,const Bool_t &i_IsoMu21_,const Bool_t &i_IsoMu21_,const Bool_t &i_IsoMu22_,const Bool_t &i_IsoMu22_,const Bool_t &i_IsoMu24_,const Bool_t &i_IsoMu27_,const Bool_t &i_IsoTkMu18_,const Bool_t &i_IsoTkMu20_,const Bool_t &i_IsoTkMu22_,const Bool_t &i_IsoTkMu22_,const Bool_t &i_IsoTkMu24_,const Bool_t &i_IsoTkMu27_,const Bool_t &i_JetE30_,const Bool_t &i_JetE30_,const Bool_t &i_JetE50_,const Bool_t &i_JetE70_,const Bool_t &i_L1SingleMu18_,const Bool_t &i_L2Mu10_,const Bool_t &i_L1SingleMuOpen_,const Bool_t &i_L1SingleMuOpen_,const Bool_t &i_L2DoubleMu23_,const Bool_t &i_L2DoubleMu28_,const Bool_t &i_L2DoubleMu38_,const Bool_t &i_L2Mu10_,const Bool_t &i_L2Mu10_,const Bool_t &i_L2Mu45_,const Bool_t &i_L2Mu40_,const Bool_t &i_LooseIsoPFTau50_,const Bool_t &i_LooseIsoPFTau50_,const Bool_t &i_LooseIsoPFTau50_,const Bool_t &i_LooseIsoPFTau50_,const Bool_t &i_LooseIsoPFTau50_,const Bool_t &i_PFTau120_,const Bool_t &i_PFTau140_,const Bool_t &i_VLooseIsoPFTau120_,const Bool_t &i_VLooseIsoPFTau140_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_Mu20_,const Bool_t &i_Mu20_,const Bool_t &i_Mu20_,const Bool_t &i_Mu20_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_Mu25_,const Bool_t &i_Mu27_,const Bool_t &i_Mu30_,const Bool_t &i_Mu30_,const Bool_t &i_Mu40_,const Bool_t &i_Mu40_,const Bool_t &i_Mu20_,const Bool_t &i_TkMu17_,const Bool_t &i_TkMu17_,const Bool_t &i_TkMu17_,const Bool_t &i_TkMu20_,const Bool_t &i_Mu24_,const Bool_t &i_TkMu24_,const Bool_t &i_Mu27_,const Bool_t &i_TkMu27_,const Bool_t &i_Mu45_,const Bool_t &i_Mu50_,const Bool_t &i_TkMu50_,const Bool_t &i_Mu38NoFiltersNoVtx_,const Bool_t &i_Mu42NoFiltersNoVtx_,const Bool_t &i_Mu28NoFiltersNoVtxDisplaced_,const Bool_t &i_Mu33NoFiltersNoVtxDisplaced_,const Bool_t &i_Mu23NoFiltersNoVtx_,const Bool_t &i_DoubleMu18NoFiltersNoVtx_,const Bool_t &i_Mu33NoFiltersNoVtxDisplaced_,const Bool_t &i_Mu33NoFiltersNoVtxDisplaced_,const Bool_t &i_Mu28NoFiltersNoVtx_,const Bool_t &i_Mu38NoFiltersNoVtxDisplaced_,const Bool_t &i_Mu38NoFiltersNoVtxDisplaced_,const Bool_t &i_Mu38NoFiltersNoVtx_,const Bool_t &i_Mu28NoFiltersNoVtx_,const Bool_t &i_PFHT300_,const Bool_t &i_PFHT300_,const Bool_t &i_PFHT550_,const Bool_t &i_PFHT650_,const Bool_t &i_PFHT750_,const Bool_t &i_PFHT750_,const Bool_t &i_PFHT750_,const Bool_t &i_PFHT800_,const Bool_t &i_PFHT850_,const Bool_t &i_PFJet15_,const Bool_t &i_PFJet25_,const Bool_t &i_DiPFJet15_,const Bool_t &i_DiPFJet25_,const Bool_t &i_DiPFJet15_,const Bool_t &i_DiPFJet25_,const Bool_t &i_DiPFJetAve15_,const Bool_t &i_DiPFJetAve25_,const Bool_t &i_DiPFJetAve35_,const Bool_t &i_AK8PFJet40_,const Bool_t &i_AK8PFJet60_,const Bool_t &i_AK8PFJet80_,const Bool_t &i_AK8PFJet140_,const Bool_t &i_AK8PFJet200_,const Bool_t &i_AK8PFJet260_,const Bool_t &i_AK8PFJet320_,const Bool_t &i_AK8PFJet400_,const Bool_t &i_AK8PFJet450_,const Bool_t &i_AK8PFJet500_,const Bool_t &i_PFJet40_,const Bool_t &i_PFJet60_,const Bool_t &i_PFJet80_,const Bool_t &i_PFJet140_,const Bool_t &i_PFJet200_,const Bool_t &i_PFJet260_,const Bool_t &i_PFJet320_,const Bool_t &i_PFJet400_,const Bool_t &i_PFJet450_,const Bool_t &i_PFJet500_,const Bool_t &i_DiPFJetAve40_,const Bool_t &i_DiPFJetAve60_,const Bool_t &i_DiPFJetAve80_,const Bool_t &i_DiPFJetAve140_,const Bool_t &i_DiPFJetAve200_,const Bool_t &i_DiPFJetAve260_,const Bool_t &i_DiPFJetAve320_,const Bool_t &i_DiPFJetAve400_,const Bool_t &i_DiPFJetAve500_,const Bool_t &i_DiPFJetAve60_,const Bool_t &i_DiPFJetAve80_,const Bool_t &i_DiPFJetAve100_,const Bool_t &i_DiPFJetAve160_,const Bool_t &i_DiPFJetAve220_,const Bool_t &i_DiPFJetAve300_,const Bool_t &i_DiPFJet40_,const Bool_t &i_DiPFJet40_,const Bool_t &i_DiCentralPFJet170_,const Bool_t &i_SingleCentralPFJet170_,const Bool_t &i_DiCentralPFJet170_,const Bool_t &i_DiCentralPFJet220_,const Bool_t &i_DiCentralPFJet330_,const Bool_t &i_DiCentralPFJet430_,const Bool_t &i_PFHT125_,const Bool_t &i_PFHT200_,const Bool_t &i_PFHT250_,const Bool_t &i_PFHT300_,const Bool_t &i_PFHT350_,const Bool_t &i_PFHT400_,const Bool_t &i_PFHT475_,const Bool_t &i_PFHT600_,const Bool_t &i_PFHT650_,const Bool_t &i_PFHT800_,const Bool_t &i_PFHT900_,const Bool_t &i_PFHT200_,const Bool_t &i_PFHT200_,const Bool_t &i_PFHT200_,const Bool_t &i_PFHT250_,const Bool_t &i_PFHT250_,const Bool_t &i_PFHT300_,const Bool_t &i_PFHT300_,const Bool_t &i_PFHT350_,const Bool_t &i_PFHT350_,const Bool_t &i_PFHT400_,const Bool_t &i_PFHT400_,const Bool_t &i_MET60_,const Bool_t &i_MET75_,const Bool_t &i_MET90_,const Bool_t &i_PFMET120_,const Bool_t &i_PFMET120_,const Bool_t &i_PFMET170_,const Bool_t &i_PFMET170_,const Bool_t &i_PFMET170_,const Bool_t &i_PFMET170_,const Bool_t &i_PFMET170_,const Bool_t &i_PFMET170_,const Bool_t &i_PFMETTypeOne190_,const Bool_t &i_PFMET90_,const Bool_t &i_PFMET100_,const Bool_t &i_PFMET100_,const Bool_t &i_PFMET110_,const Bool_t &i_PFMET120_,const Bool_t &i_CaloMHTNoPU90_,const Bool_t &i_CaloMHTNoPU90_,const Bool_t &i_QuadPFJet_,const Bool_t &i_QuadPFJet_,const Bool_t &i_QuadPFJet_,const Bool_t &i_QuadPFJet_,const Bool_t &i_QuadPFJet_,const Bool_t &i_L1_,const Bool_t &i_QuadJet45_,const Bool_t &i_QuadJet45_,const Bool_t &i_DoubleJet90_,const Bool_t &i_DoubleJet90_,const Bool_t &i_DoubleJetsC100_,const Bool_t &i_DoubleJetsC100_,const Bool_t &i_DoubleJetsC112_,const Bool_t &i_DoubleJetsC112_,const Bool_t &i_DoubleJetsC100_,const Bool_t &i_DoubleJetsC100_,const Bool_t &i_DoubleJetsC100_,const Bool_t &i_DoubleJetsC100_,const Bool_t &i_Photon135_,const Bool_t &i_Photon20_,const Bool_t &i_Photon22_,const Bool_t &i_Photon22_,const Bool_t &i_Photon250_,const Bool_t &i_Photon300_,const Bool_t &i_Photon26_,const Bool_t &i_Photon36_,const Bool_t &i_Photon36_,const Bool_t &i_Photon36_,const Bool_t &i_Photon50_,const Bool_t &i_Photon50_,const Bool_t &i_Photon75_,const Bool_t &i_Photon75_,const Bool_t &i_Photon90_,const Bool_t &i_Photon90_,const Bool_t &i_Photon120_,const Bool_t &i_Photon120_,const Bool_t &i_Mu8_,const Bool_t &i_Mu17_,const Bool_t &i_Ele8_,const Bool_t &i_Ele12_,const Bool_t &i_Ele17_,const Bool_t &i_Ele23_,const Bool_t &i_BTagMu_,const Bool_t &i_BTagMu_,const Bool_t &i_BTagMu_,const Bool_t &i_BTagMu_,const Bool_t &i_BTagMu_,const Bool_t &i_BTagMu_,const Bool_t &i_BTagMu_,const Bool_t &i_Ele23_,const Bool_t &i_Ele23_,const Bool_t &i_Ele17_,const Bool_t &i_Ele16_,const Bool_t &i_Mu8_,const Bool_t &i_Mu8_,const Bool_t &i_Mu8_,const Bool_t &i_Mu12_,const Bool_t &i_Mu12_,const Bool_t &i_Mu17_,const Bool_t &i_Mu23_,const Bool_t &i_Mu23_,const Bool_t &i_Mu23_,const Bool_t &i_Mu23_,const Bool_t &i_Mu30_,const Bool_t &i_Mu33_,const Bool_t &i_Mu37_,const Bool_t &i_Mu27_,const Bool_t &i_Mu8_,const Bool_t &i_Mu12_,const Bool_t &i_Mu12_,const Bool_t &i_Mu12_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_Mu17_,const Bool_t &i_DiMu9_,const Bool_t &i_TripleMu_,const Bool_t &i_TripleMu_,const Bool_t &i_Mu3er_,const Bool_t &i_Mu6_,const Bool_t &i_Mu6_,const Bool_t &i_Mu14er_,const Bool_t &i_Ele17_,const Bool_t &i_Ele23_,const Bool_t &i_Ele12_,const Bool_t &i_Ele17_,const Bool_t &i_Ele17_,const Bool_t &i_Ele23_,const Bool_t &i_PFHT650_,const Bool_t &i_PFHT650_,const Bool_t &i_Photon22_,const Bool_t &i_Photon30_,const Bool_t &i_Photon36_,const Bool_t &i_Photon50_,const Bool_t &i_Photon75_,const Bool_t &i_Photon90_,const Bool_t &i_Photon120_,const Bool_t &i_Photon175_,const Bool_t &i_Photon165_,const Bool_t &i_Photon22_,const Bool_t &i_Photon30_,const Bool_t &i_Photon36_,const Bool_t &i_Photon50_,const Bool_t &i_Photon75_,const Bool_t &i_Photon90_,const Bool_t &i_Photon120_,const Bool_t &i_Photon165_,const Bool_t &i_Diphoton30_,const Bool_t &i_Diphoton30_,const Bool_t &i_Diphoton30PV_,const Bool_t &i_Diphoton30_,const Bool_t &i_Diphoton30EB_,const Bool_t &i_Dimuon0_,const Bool_t &i_Dimuon0_,const Bool_t &i_QuadMuon0_,const Bool_t &i_QuadMuon0_,const Bool_t &i_Rsq0p25_,const Bool_t &i_RsqMR240_,const Bool_t &i_RsqMR240_,const Bool_t &i_Rsq0p25_,const Bool_t &i_Rsq0p30_,const Bool_t &i_RsqMR240_,const Bool_t &i_RsqMR240_,const Bool_t &i_RsqMR270_,const Bool_t &i_RsqMR270_,const Bool_t &i_Rsq0p02_,const Bool_t &i_Rsq0p02_,const Bool_t &i_Rsq0p02_,const Bool_t &i_Rsq0p02_,const Bool_t &i_Rsq0p02_,const Bool_t &i_HT200_,const Bool_t &i_HT250_,const Bool_t &i_HT350_,const Bool_t &i_HT350_,const Bool_t &i_HT350_,const Bool_t &i_HT350_,const Bool_t &i_HT400_,const Bool_t &i_HT500_,const Bool_t &i_HT550_,const Bool_t &i_HT550_,const Bool_t &i_HT650_,const Bool_t &i_HT750_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_VBF_,const Bool_t &i_PFMETNoMu90_,const Bool_t &i_PFMETNoMu100_,const Bool_t &i_PFMETNoMu110_,const Bool_t &i_PFMETNoMu120_,const Bool_t &i_MonoCentralPFJet80_,const Bool_t &i_MonoCentralPFJet80_,const Bool_t &i_MonoCentralPFJet80_,const Bool_t &i_MonoCentralPFJet80_,const Bool_t &i_Ele27_,const Bool_t &i_Photon90_,const Bool_t &i_DoubleMu8_,const Bool_t &i_Mu8_,const Bool_t &i_DoubleEle8_,const Bool_t &i_DoubleMu8_,const Bool_t &i_Mu8_,const Bool_t &i_DoubleEle8_,const Bool_t &i_Mu10_,const Bool_t &i_DoubleMu3_,const Bool_t &i_Ele10_,const Bool_t &i_Ele15_,const Bool_t &i_Ele15_,const Bool_t &i_Ele15_,const Bool_t &i_Ele15_,const Bool_t &i_Ele15_,const Bool_t &i_Ele15_,const Bool_t &i_Ele50_,const Bool_t &i_Mu8_,const Bool_t &i_Mu10_,const Bool_t &i_Mu15_,const Bool_t &i_Mu15_,const Bool_t &i_Mu15_,const Bool_t &i_Mu15_,const Bool_t &i_Mu15_,const Bool_t &i_Mu15_,const Bool_t &i_Mu50_,const Bool_t &i_Dimuon16_,const Bool_t &i_Dimuon10_,const Bool_t &i_Dimuon8_,const Bool_t &i_Dimuon8_,const Bool_t &i_Dimuon0_,const Bool_t &i_Mu16_,const Bool_t &i_Mu16_,const Bool_t &i_TrkMu15_,const Bool_t &i_TrkMu17_,const Bool_t &i_Mu8_,const Bool_t &i_Mu17_,const Bool_t &i_Mu3_,const Bool_t &i_Ele8_,const Bool_t &i_Ele12_,const Bool_t &i_Ele17_,const Bool_t &i_Ele23_,const Bool_t &i_Ele50_,const Bool_t &i_Ele50_,const Bool_t &i_PFHT400_,const Bool_t &i_PFHT450_,const Bool_t &i_PFHT400_,const Bool_t &i_PFHT450_,const Bool_t &i_Ele115_,const Bool_t &i_Mu55_,const Bool_t &i_Photon42_,const Bool_t &i_Photon90_,const Bool_t &i_PixelTracks_,const Bool_t &i_PixelTracks_,const Bool_t &i_PixelTracks_,const Bool_t &i_PixelTracks_,const Bool_t &i_PixelTracks_,const Bool_t &i_FullTracks_,const Bool_t &i_FullTracks_,const Bool_t &i_FullTracks_,const Bool_t &i_FullTracks_,const Bool_t &i_ECALHT800_,const Bool_t &i_DiSC30_,const Bool_t &i_Photon125_,const Bool_t &i_MET100_,const Bool_t &i_MET150_,const Bool_t &i_MET200_,const Bool_t &i_Ele27_,const Bool_t &i_L1FatEvents_,const Bool_t &i_Physics_,const Bool_t &i_L1FatEvents_,const Bool_t &i_L1FatEvents_,const Bool_t &i_L1FatEvents_,const Bool_t &i_L1FatEvents_,const Bool_t &i_Random_,const Bool_t &i_ZeroBias_,const Bool_t &i_AK4CaloJet30_,const Bool_t &i_AK4CaloJet40_,const Bool_t &i_AK4CaloJet50_,const Bool_t &i_AK4CaloJet80_,const Bool_t &i_AK4CaloJet100_,const Bool_t &i_AK4PFJet30_,const Bool_t &i_AK4PFJet50_,const Bool_t &i_AK4PFJet80_,const Bool_t &i_AK4PFJet100_,const Bool_t &i_HISinglePhoton10_,const Bool_t &i_HISinglePhoton15_,const Bool_t &i_HISinglePhoton20_,const Bool_t &i_HISinglePhoton40_,const Bool_t &i_HISinglePhoton60_,const Bool_t &i_EcalCalibration_,const Bool_t &i_HcalCalibration_,const Bool_t &i_GlobalRunHPDNoise_,const Bool_t &i_L1BptxMinus_,const Bool_t &i_L1BptxPlus_,const Bool_t &i_L1NotBptxOR_,const Bool_t &i_L1BeamGasMinus_,const Bool_t &i_L1BeamGasPlus_,const Bool_t &i_L1BptxXOR_,const Bool_t &i_L1MinimumBiasHF_,const Bool_t &i_L1MinimumBiasHF_,const Bool_t &i_HcalNZS_,const Bool_t &i_HcalPhiSym_,const Bool_t &i_HcalIsolatedbunch_,const Bool_t &i_ZeroBias_,const Bool_t &i_ZeroBias_,const Bool_t &i_ZeroBias_,const Bool_t &i_ZeroBias_,const Bool_t &i_ZeroBias_,const Bool_t &i_ZeroBias_,const Bool_t &i_Photon500_,const Bool_t &i_Photon600_,const Bool_t &i_Mu300_,const Bool_t &i_Mu350_,const Bool_t &i_MET250_,const Bool_t &i_MET300_,const Bool_t &i_MET600_,const Bool_t &i_MET700_,const Bool_t &i_PFMET300_,const Bool_t &i_PFMET400_,const Bool_t &i_PFMET500_,const Bool_t &i_PFMET600_,const Bool_t &i_Ele250_,const Bool_t &i_Ele300_,const Bool_t &i_HT2000_,const Bool_t &i_HT2500_,const Bool_t &i_IsoTrackHE_,const Bool_t &i_IsoTrackHB_,const Bool_t &i_HLTriggerFinalPath_):
//    
//  {}
  Hlt():
    HLTriggerFirstPath_(0),
    AK8PFJet360_(0),
    AK8PFJet400_(0),
    AK8PFHT750_(0),
    AK8PFHT800_(0),
    AK8DiPFJet300_(0),
    AK8DiPFJet280_(0),
    AK8DiPFJet300_(0),
    AK8DiPFJet300_(0),
    AK8PFHT700_(0),
    AK8PFHT650_(0),
    AK8PFHT600_(0),
    AK8DiPFJet280_(0),
    AK8DiPFJet250_(0),
    AK8DiPFJet280_(0),
    AK8DiPFJet250_(0),
    CaloJet260_(0),
    CaloJet500_(0),
    Dimuon13_(0),
    Dimuon13_(0),
    Dimuon20_(0),
    DoubleEle24_(0),
    DoubleEle25_(0),
    DoubleEle33_(0),
    DoubleEle33_(0),
    DoubleEle33_(0),
    DoubleEle33_(0),
    DoubleMediumCombinedIsoPFTau35_(0),
    DoubleTightCombinedIsoPFTau35_(0),
    DoubleMediumCombinedIsoPFTau40_(0),
    DoubleTightCombinedIsoPFTau40_(0),
    DoubleMediumCombinedIsoPFTau40_(0),
    DoubleTightCombinedIsoPFTau40_(0),
    DoubleMediumIsoPFTau35_(0),
    DoubleMediumIsoPFTau40_(0),
    DoubleMediumIsoPFTau40_(0),
    DoubleEle37_(0),
    DoubleMu33NoFiltersNoVtx_(0),
    DoubleMu38NoFiltersNoVtx_(0),
    DoubleMu23NoFiltersNoVtxDisplaced_(0),
    DoubleMu28NoFiltersNoVtxDisplaced_(0),
    DoubleMu0_(0),
    DoubleMu4_(0),
    DoubleMu4_(0),
    DoubleMu4_(0),
    DoubleMu4_(0),
    DoubleMu3_(0),
    DoubleMu4_(0),
    Mu7p5_(0),
    Mu7p5_(0),
    Mu7p5_(0),
    Mu7p5_(0),
    Mu7p5_(0),
    Mu7p5_(0),
    Mu7p5_(0),
    Mu7p5_(0),
    Dimuon0er16_(0),
    Dimuon0er16_(0),
    Dimuon6_(0),
    Photon150_(0),
    Photon90_(0),
    HT250_(0),
    DoublePhoton60_(0),
    DoublePhoton85_(0),
    Ele17_(0),
    Ele20_(0),
    Ele22_(0),
    Ele22_(0),
    Ele22_(0),
    Ele23_(0),
    Ele23_(0),
    Ele24_(0),
    Ele24_(0),
    Ele24_(0),
    Ele24_(0),
    Ele25_(0),
    Ele25_(0),
    Ele25_(0),
    Ele27_(0),
    Ele27_(0),
    Ele27_(0),
    Ele27_(0),
    Ele27_(0),
    Ele27_(0),
    Ele27_(0),
    Ele30_(0),
    Ele30_(0),
    Ele30_(0),
    Ele32_(0),
    Ele32_(0),
    Ele32_(0),
    Ele32_(0),
    Ele35_(0),
    Ele35_(0),
    Ele36_(0),
    Ele45_(0),
    Ele45_(0),
    Ele45_(0),
    Ele105_(0),
    Ele30WP60_(0),
    Ele30WP60_(0),
    HT200_(0),
    HT275_(0),
    HT325_(0),
    HT425_(0),
    HT575_(0),
    HT410to430_(0),
    HT430to450_(0),
    HT450to470_(0),
    HT470to500_(0),
    HT500to550_(0),
    HT550to650_(0),
    HT650_(0),
    Mu16_(0),
    IsoMu16_(0),
    IsoMu16_(0),
    IsoMu17_(0),
    IsoMu17_(0),
    IsoMu17_(0),
    DoubleIsoMu17_(0),
    DoubleIsoMu17_(0),
    IsoMu18_(0),
    IsoMu19_(0),
    IsoMu19_(0),
    IsoMu19_(0),
    IsoMu19_(0),
    IsoMu19_(0),
    IsoMu19_(0),
    IsoMu21_(0),
    IsoMu21_(0),
    IsoMu20_(0),
    IsoMu21_(0),
    IsoMu21_(0),
    IsoMu21_(0),
    IsoMu22_(0),
    IsoMu22_(0),
    IsoMu24_(0),
    IsoMu27_(0),
    IsoTkMu18_(0),
    IsoTkMu20_(0),
    IsoTkMu22_(0),
    IsoTkMu22_(0),
    IsoTkMu24_(0),
    IsoTkMu27_(0),
    JetE30_(0),
    JetE30_(0),
    JetE50_(0),
    JetE70_(0),
    L1SingleMu18_(0),
    L2Mu10_(0),
    L1SingleMuOpen_(0),
    L1SingleMuOpen_(0),
    L2DoubleMu23_(0),
    L2DoubleMu28_(0),
    L2DoubleMu38_(0),
    L2Mu10_(0),
    L2Mu10_(0),
    L2Mu45_(0),
    L2Mu40_(0),
    LooseIsoPFTau50_(0),
    LooseIsoPFTau50_(0),
    LooseIsoPFTau50_(0),
    LooseIsoPFTau50_(0),
    LooseIsoPFTau50_(0),
    PFTau120_(0),
    PFTau140_(0),
    VLooseIsoPFTau120_(0),
    VLooseIsoPFTau140_(0),
    Mu17_(0),
    Mu17_(0),
    Mu17_(0),
    Mu17_(0),
    Mu20_(0),
    Mu20_(0),
    Mu20_(0),
    Mu20_(0),
    Mu17_(0),
    Mu17_(0),
    Mu17_(0),
    Mu17_(0),
    Mu17_(0),
    Mu25_(0),
    Mu27_(0),
    Mu30_(0),
    Mu30_(0),
    Mu40_(0),
    Mu40_(0),
    Mu20_(0),
    TkMu17_(0),
    TkMu17_(0),
    TkMu17_(0),
    TkMu20_(0),
    Mu24_(0),
    TkMu24_(0),
    Mu27_(0),
    TkMu27_(0),
    Mu45_(0),
    Mu50_(0),
    TkMu50_(0),
    Mu38NoFiltersNoVtx_(0),
    Mu42NoFiltersNoVtx_(0),
    Mu28NoFiltersNoVtxDisplaced_(0),
    Mu33NoFiltersNoVtxDisplaced_(0),
    Mu23NoFiltersNoVtx_(0),
    DoubleMu18NoFiltersNoVtx_(0),
    Mu33NoFiltersNoVtxDisplaced_(0),
    Mu33NoFiltersNoVtxDisplaced_(0),
    Mu28NoFiltersNoVtx_(0),
    Mu38NoFiltersNoVtxDisplaced_(0),
    Mu38NoFiltersNoVtxDisplaced_(0),
    Mu38NoFiltersNoVtx_(0),
    Mu28NoFiltersNoVtx_(0),
    PFHT300_(0),
    PFHT300_(0),
    PFHT550_(0),
    PFHT650_(0),
    PFHT750_(0),
    PFHT750_(0),
    PFHT750_(0),
    PFHT800_(0),
    PFHT850_(0),
    PFJet15_(0),
    PFJet25_(0),
    DiPFJet15_(0),
    DiPFJet25_(0),
    DiPFJet15_(0),
    DiPFJet25_(0),
    DiPFJetAve15_(0),
    DiPFJetAve25_(0),
    DiPFJetAve35_(0),
    AK8PFJet40_(0),
    AK8PFJet60_(0),
    AK8PFJet80_(0),
    AK8PFJet140_(0),
    AK8PFJet200_(0),
    AK8PFJet260_(0),
    AK8PFJet320_(0),
    AK8PFJet400_(0),
    AK8PFJet450_(0),
    AK8PFJet500_(0),
    PFJet40_(0),
    PFJet60_(0),
    PFJet80_(0),
    PFJet140_(0),
    PFJet200_(0),
    PFJet260_(0),
    PFJet320_(0),
    PFJet400_(0),
    PFJet450_(0),
    PFJet500_(0),
    DiPFJetAve40_(0),
    DiPFJetAve60_(0),
    DiPFJetAve80_(0),
    DiPFJetAve140_(0),
    DiPFJetAve200_(0),
    DiPFJetAve260_(0),
    DiPFJetAve320_(0),
    DiPFJetAve400_(0),
    DiPFJetAve500_(0),
    DiPFJetAve60_(0),
    DiPFJetAve80_(0),
    DiPFJetAve100_(0),
    DiPFJetAve160_(0),
    DiPFJetAve220_(0),
    DiPFJetAve300_(0),
    DiPFJet40_(0),
    DiPFJet40_(0),
    DiCentralPFJet170_(0),
    SingleCentralPFJet170_(0),
    DiCentralPFJet170_(0),
    DiCentralPFJet220_(0),
    DiCentralPFJet330_(0),
    DiCentralPFJet430_(0),
    PFHT125_(0),
    PFHT200_(0),
    PFHT250_(0),
    PFHT300_(0),
    PFHT350_(0),
    PFHT400_(0),
    PFHT475_(0),
    PFHT600_(0),
    PFHT650_(0),
    PFHT800_(0),
    PFHT900_(0),
    PFHT200_(0),
    PFHT200_(0),
    PFHT200_(0),
    PFHT250_(0),
    PFHT250_(0),
    PFHT300_(0),
    PFHT300_(0),
    PFHT350_(0),
    PFHT350_(0),
    PFHT400_(0),
    PFHT400_(0),
    MET60_(0),
    MET75_(0),
    MET90_(0),
    PFMET120_(0),
    PFMET120_(0),
    PFMET170_(0),
    PFMET170_(0),
    PFMET170_(0),
    PFMET170_(0),
    PFMET170_(0),
    PFMET170_(0),
    PFMETTypeOne190_(0),
    PFMET90_(0),
    PFMET100_(0),
    PFMET100_(0),
    PFMET110_(0),
    PFMET120_(0),
    CaloMHTNoPU90_(0),
    CaloMHTNoPU90_(0),
    QuadPFJet_(0),
    QuadPFJet_(0),
    QuadPFJet_(0),
    QuadPFJet_(0),
    QuadPFJet_(0),
    L1_(0),
    QuadJet45_(0),
    QuadJet45_(0),
    DoubleJet90_(0),
    DoubleJet90_(0),
    DoubleJetsC100_(0),
    DoubleJetsC100_(0),
    DoubleJetsC112_(0),
    DoubleJetsC112_(0),
    DoubleJetsC100_(0),
    DoubleJetsC100_(0),
    DoubleJetsC100_(0),
    DoubleJetsC100_(0),
    Photon135_(0),
    Photon20_(0),
    Photon22_(0),
    Photon22_(0),
    Photon250_(0),
    Photon300_(0),
    Photon26_(0),
    Photon36_(0),
    Photon36_(0),
    Photon36_(0),
    Photon50_(0),
    Photon50_(0),
    Photon75_(0),
    Photon75_(0),
    Photon90_(0),
    Photon90_(0),
    Photon120_(0),
    Photon120_(0),
    Mu8_(0),
    Mu17_(0),
    Ele8_(0),
    Ele12_(0),
    Ele17_(0),
    Ele23_(0),
    BTagMu_(0),
    BTagMu_(0),
    BTagMu_(0),
    BTagMu_(0),
    BTagMu_(0),
    BTagMu_(0),
    BTagMu_(0),
    Ele23_(0),
    Ele23_(0),
    Ele17_(0),
    Ele16_(0),
    Mu8_(0),
    Mu8_(0),
    Mu8_(0),
    Mu12_(0),
    Mu12_(0),
    Mu17_(0),
    Mu23_(0),
    Mu23_(0),
    Mu23_(0),
    Mu23_(0),
    Mu30_(0),
    Mu33_(0),
    Mu37_(0),
    Mu27_(0),
    Mu8_(0),
    Mu12_(0),
    Mu12_(0),
    Mu12_(0),
    Mu17_(0),
    Mu17_(0),
    Mu17_(0),
    DiMu9_(0),
    TripleMu_(0),
    TripleMu_(0),
    Mu3er_(0),
    Mu6_(0),
    Mu6_(0),
    Mu14er_(0),
    Ele17_(0),
    Ele23_(0),
    Ele12_(0),
    Ele17_(0),
    Ele17_(0),
    Ele23_(0),
    PFHT650_(0),
    PFHT650_(0),
    Photon22_(0),
    Photon30_(0),
    Photon36_(0),
    Photon50_(0),
    Photon75_(0),
    Photon90_(0),
    Photon120_(0),
    Photon175_(0),
    Photon165_(0),
    Photon22_(0),
    Photon30_(0),
    Photon36_(0),
    Photon50_(0),
    Photon75_(0),
    Photon90_(0),
    Photon120_(0),
    Photon165_(0),
    Diphoton30_(0),
    Diphoton30_(0),
    Diphoton30PV_(0),
    Diphoton30_(0),
    Diphoton30EB_(0),
    Dimuon0_(0),
    Dimuon0_(0),
    QuadMuon0_(0),
    QuadMuon0_(0),
    Rsq0p25_(0),
    RsqMR240_(0),
    RsqMR240_(0),
    Rsq0p25_(0),
    Rsq0p30_(0),
    RsqMR240_(0),
    RsqMR240_(0),
    RsqMR270_(0),
    RsqMR270_(0),
    Rsq0p02_(0),
    Rsq0p02_(0),
    Rsq0p02_(0),
    Rsq0p02_(0),
    Rsq0p02_(0),
    HT200_(0),
    HT250_(0),
    HT350_(0),
    HT350_(0),
    HT350_(0),
    HT350_(0),
    HT400_(0),
    HT500_(0),
    HT550_(0),
    HT550_(0),
    HT650_(0),
    HT750_(0),
    VBF_(0),
    VBF_(0),
    VBF_(0),
    VBF_(0),
    VBF_(0),
    VBF_(0),
    VBF_(0),
    VBF_(0),
    VBF_(0),
    VBF_(0),
    PFMETNoMu90_(0),
    PFMETNoMu100_(0),
    PFMETNoMu110_(0),
    PFMETNoMu120_(0),
    MonoCentralPFJet80_(0),
    MonoCentralPFJet80_(0),
    MonoCentralPFJet80_(0),
    MonoCentralPFJet80_(0),
    Ele27_(0),
    Photon90_(0),
    DoubleMu8_(0),
    Mu8_(0),
    DoubleEle8_(0),
    DoubleMu8_(0),
    Mu8_(0),
    DoubleEle8_(0),
    Mu10_(0),
    DoubleMu3_(0),
    Ele10_(0),
    Ele15_(0),
    Ele15_(0),
    Ele15_(0),
    Ele15_(0),
    Ele15_(0),
    Ele15_(0),
    Ele50_(0),
    Mu8_(0),
    Mu10_(0),
    Mu15_(0),
    Mu15_(0),
    Mu15_(0),
    Mu15_(0),
    Mu15_(0),
    Mu15_(0),
    Mu50_(0),
    Dimuon16_(0),
    Dimuon10_(0),
    Dimuon8_(0),
    Dimuon8_(0),
    Dimuon0_(0),
    Mu16_(0),
    Mu16_(0),
    TrkMu15_(0),
    TrkMu17_(0),
    Mu8_(0),
    Mu17_(0),
    Mu3_(0),
    Ele8_(0),
    Ele12_(0),
    Ele17_(0),
    Ele23_(0),
    Ele50_(0),
    Ele50_(0),
    PFHT400_(0),
    PFHT450_(0),
    PFHT400_(0),
    PFHT450_(0),
    Ele115_(0),
    Mu55_(0),
    Photon42_(0),
    Photon90_(0),
    PixelTracks_(0),
    PixelTracks_(0),
    PixelTracks_(0),
    PixelTracks_(0),
    PixelTracks_(0),
    FullTracks_(0),
    FullTracks_(0),
    FullTracks_(0),
    FullTracks_(0),
    ECALHT800_(0),
    DiSC30_(0),
    Photon125_(0),
    MET100_(0),
    MET150_(0),
    MET200_(0),
    Ele27_(0),
    L1FatEvents_(0),
    Physics_(0),
    L1FatEvents_(0),
    L1FatEvents_(0),
    L1FatEvents_(0),
    L1FatEvents_(0),
    Random_(0),
    ZeroBias_(0),
    AK4CaloJet30_(0),
    AK4CaloJet40_(0),
    AK4CaloJet50_(0),
    AK4CaloJet80_(0),
    AK4CaloJet100_(0),
    AK4PFJet30_(0),
    AK4PFJet50_(0),
    AK4PFJet80_(0),
    AK4PFJet100_(0),
    HISinglePhoton10_(0),
    HISinglePhoton15_(0),
    HISinglePhoton20_(0),
    HISinglePhoton40_(0),
    HISinglePhoton60_(0),
    EcalCalibration_(0),
    HcalCalibration_(0),
    GlobalRunHPDNoise_(0),
    L1BptxMinus_(0),
    L1BptxPlus_(0),
    L1NotBptxOR_(0),
    L1BeamGasMinus_(0),
    L1BeamGasPlus_(0),
    L1BptxXOR_(0),
    L1MinimumBiasHF_(0),
    L1MinimumBiasHF_(0),
    HcalNZS_(0),
    HcalPhiSym_(0),
    HcalIsolatedbunch_(0),
    ZeroBias_(0),
    ZeroBias_(0),
    ZeroBias_(0),
    ZeroBias_(0),
    ZeroBias_(0),
    ZeroBias_(0),
    Photon500_(0),
    Photon600_(0),
    Mu300_(0),
    Mu350_(0),
    MET250_(0),
    MET300_(0),
    MET600_(0),
    MET700_(0),
    PFMET300_(0),
    PFMET400_(0),
    PFMET500_(0),
    PFMET600_(0),
    Ele250_(0),
    Ele300_(0),
    HT2000_(0),
    HT2500_(0),
    IsoTrackHE_(0),
    IsoTrackHB_(0),
    HLTriggerFinalPath_(0)
  {}
  Bool_t HLTriggerFirstPath() const {return HLTriggerFirstPath_;}
  Bool_t AK8PFJet360() const {return AK8PFJet360_;}
  Bool_t AK8PFJet400() const {return AK8PFJet400_;}
  Bool_t AK8PFHT750() const {return AK8PFHT750_;}
  Bool_t AK8PFHT800() const {return AK8PFHT800_;}
  Bool_t AK8DiPFJet300() const {return AK8DiPFJet300_;}
  Bool_t AK8DiPFJet280() const {return AK8DiPFJet280_;}
  Bool_t AK8DiPFJet300() const {return AK8DiPFJet300_;}
  Bool_t AK8DiPFJet300() const {return AK8DiPFJet300_;}
  Bool_t AK8PFHT700() const {return AK8PFHT700_;}
  Bool_t AK8PFHT650() const {return AK8PFHT650_;}
  Bool_t AK8PFHT600() const {return AK8PFHT600_;}
  Bool_t AK8DiPFJet280() const {return AK8DiPFJet280_;}
  Bool_t AK8DiPFJet250() const {return AK8DiPFJet250_;}
  Bool_t AK8DiPFJet280() const {return AK8DiPFJet280_;}
  Bool_t AK8DiPFJet250() const {return AK8DiPFJet250_;}
  Bool_t CaloJet260() const {return CaloJet260_;}
  Bool_t CaloJet500() const {return CaloJet500_;}
  Bool_t Dimuon13() const {return Dimuon13_;}
  Bool_t Dimuon13() const {return Dimuon13_;}
  Bool_t Dimuon20() const {return Dimuon20_;}
  Bool_t DoubleEle24() const {return DoubleEle24_;}
  Bool_t DoubleEle25() const {return DoubleEle25_;}
  Bool_t DoubleEle33() const {return DoubleEle33_;}
  Bool_t DoubleEle33() const {return DoubleEle33_;}
  Bool_t DoubleEle33() const {return DoubleEle33_;}
  Bool_t DoubleEle33() const {return DoubleEle33_;}
  Bool_t DoubleMediumCombinedIsoPFTau35() const {return DoubleMediumCombinedIsoPFTau35_;}
  Bool_t DoubleTightCombinedIsoPFTau35() const {return DoubleTightCombinedIsoPFTau35_;}
  Bool_t DoubleMediumCombinedIsoPFTau40() const {return DoubleMediumCombinedIsoPFTau40_;}
  Bool_t DoubleTightCombinedIsoPFTau40() const {return DoubleTightCombinedIsoPFTau40_;}
  Bool_t DoubleMediumCombinedIsoPFTau40() const {return DoubleMediumCombinedIsoPFTau40_;}
  Bool_t DoubleTightCombinedIsoPFTau40() const {return DoubleTightCombinedIsoPFTau40_;}
  Bool_t DoubleMediumIsoPFTau35() const {return DoubleMediumIsoPFTau35_;}
  Bool_t DoubleMediumIsoPFTau40() const {return DoubleMediumIsoPFTau40_;}
  Bool_t DoubleMediumIsoPFTau40() const {return DoubleMediumIsoPFTau40_;}
  Bool_t DoubleEle37() const {return DoubleEle37_;}
  Bool_t DoubleMu33NoFiltersNoVtx() const {return DoubleMu33NoFiltersNoVtx_;}
  Bool_t DoubleMu38NoFiltersNoVtx() const {return DoubleMu38NoFiltersNoVtx_;}
  Bool_t DoubleMu23NoFiltersNoVtxDisplaced() const {return DoubleMu23NoFiltersNoVtxDisplaced_;}
  Bool_t DoubleMu28NoFiltersNoVtxDisplaced() const {return DoubleMu28NoFiltersNoVtxDisplaced_;}
  Bool_t DoubleMu0() const {return DoubleMu0_;}
  Bool_t DoubleMu4() const {return DoubleMu4_;}
  Bool_t DoubleMu4() const {return DoubleMu4_;}
  Bool_t DoubleMu4() const {return DoubleMu4_;}
  Bool_t DoubleMu4() const {return DoubleMu4_;}
  Bool_t DoubleMu3() const {return DoubleMu3_;}
  Bool_t DoubleMu4() const {return DoubleMu4_;}
  Bool_t Mu7p5() const {return Mu7p5_;}
  Bool_t Mu7p5() const {return Mu7p5_;}
  Bool_t Mu7p5() const {return Mu7p5_;}
  Bool_t Mu7p5() const {return Mu7p5_;}
  Bool_t Mu7p5() const {return Mu7p5_;}
  Bool_t Mu7p5() const {return Mu7p5_;}
  Bool_t Mu7p5() const {return Mu7p5_;}
  Bool_t Mu7p5() const {return Mu7p5_;}
  Bool_t Dimuon0er16() const {return Dimuon0er16_;}
  Bool_t Dimuon0er16() const {return Dimuon0er16_;}
  Bool_t Dimuon6() const {return Dimuon6_;}
  Bool_t Photon150() const {return Photon150_;}
  Bool_t Photon90() const {return Photon90_;}
  Bool_t HT250() const {return HT250_;}
  Bool_t DoublePhoton60() const {return DoublePhoton60_;}
  Bool_t DoublePhoton85() const {return DoublePhoton85_;}
  Bool_t Ele17() const {return Ele17_;}
  Bool_t Ele20() const {return Ele20_;}
  Bool_t Ele22() const {return Ele22_;}
  Bool_t Ele22() const {return Ele22_;}
  Bool_t Ele22() const {return Ele22_;}
  Bool_t Ele23() const {return Ele23_;}
  Bool_t Ele23() const {return Ele23_;}
  Bool_t Ele24() const {return Ele24_;}
  Bool_t Ele24() const {return Ele24_;}
  Bool_t Ele24() const {return Ele24_;}
  Bool_t Ele24() const {return Ele24_;}
  Bool_t Ele25() const {return Ele25_;}
  Bool_t Ele25() const {return Ele25_;}
  Bool_t Ele25() const {return Ele25_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t Ele30() const {return Ele30_;}
  Bool_t Ele30() const {return Ele30_;}
  Bool_t Ele30() const {return Ele30_;}
  Bool_t Ele32() const {return Ele32_;}
  Bool_t Ele32() const {return Ele32_;}
  Bool_t Ele32() const {return Ele32_;}
  Bool_t Ele32() const {return Ele32_;}
  Bool_t Ele35() const {return Ele35_;}
  Bool_t Ele35() const {return Ele35_;}
  Bool_t Ele36() const {return Ele36_;}
  Bool_t Ele45() const {return Ele45_;}
  Bool_t Ele45() const {return Ele45_;}
  Bool_t Ele45() const {return Ele45_;}
  Bool_t Ele105() const {return Ele105_;}
  Bool_t Ele30WP60() const {return Ele30WP60_;}
  Bool_t Ele30WP60() const {return Ele30WP60_;}
  Bool_t HT200() const {return HT200_;}
  Bool_t HT275() const {return HT275_;}
  Bool_t HT325() const {return HT325_;}
  Bool_t HT425() const {return HT425_;}
  Bool_t HT575() const {return HT575_;}
  Bool_t HT410to430() const {return HT410to430_;}
  Bool_t HT430to450() const {return HT430to450_;}
  Bool_t HT450to470() const {return HT450to470_;}
  Bool_t HT470to500() const {return HT470to500_;}
  Bool_t HT500to550() const {return HT500to550_;}
  Bool_t HT550to650() const {return HT550to650_;}
  Bool_t HT650() const {return HT650_;}
  Bool_t Mu16() const {return Mu16_;}
  Bool_t IsoMu16() const {return IsoMu16_;}
  Bool_t IsoMu16() const {return IsoMu16_;}
  Bool_t IsoMu17() const {return IsoMu17_;}
  Bool_t IsoMu17() const {return IsoMu17_;}
  Bool_t IsoMu17() const {return IsoMu17_;}
  Bool_t DoubleIsoMu17() const {return DoubleIsoMu17_;}
  Bool_t DoubleIsoMu17() const {return DoubleIsoMu17_;}
  Bool_t IsoMu18() const {return IsoMu18_;}
  Bool_t IsoMu19() const {return IsoMu19_;}
  Bool_t IsoMu19() const {return IsoMu19_;}
  Bool_t IsoMu19() const {return IsoMu19_;}
  Bool_t IsoMu19() const {return IsoMu19_;}
  Bool_t IsoMu19() const {return IsoMu19_;}
  Bool_t IsoMu19() const {return IsoMu19_;}
  Bool_t IsoMu21() const {return IsoMu21_;}
  Bool_t IsoMu21() const {return IsoMu21_;}
  Bool_t IsoMu20() const {return IsoMu20_;}
  Bool_t IsoMu21() const {return IsoMu21_;}
  Bool_t IsoMu21() const {return IsoMu21_;}
  Bool_t IsoMu21() const {return IsoMu21_;}
  Bool_t IsoMu22() const {return IsoMu22_;}
  Bool_t IsoMu22() const {return IsoMu22_;}
  Bool_t IsoMu24() const {return IsoMu24_;}
  Bool_t IsoMu27() const {return IsoMu27_;}
  Bool_t IsoTkMu18() const {return IsoTkMu18_;}
  Bool_t IsoTkMu20() const {return IsoTkMu20_;}
  Bool_t IsoTkMu22() const {return IsoTkMu22_;}
  Bool_t IsoTkMu22() const {return IsoTkMu22_;}
  Bool_t IsoTkMu24() const {return IsoTkMu24_;}
  Bool_t IsoTkMu27() const {return IsoTkMu27_;}
  Bool_t JetE30() const {return JetE30_;}
  Bool_t JetE30() const {return JetE30_;}
  Bool_t JetE50() const {return JetE50_;}
  Bool_t JetE70() const {return JetE70_;}
  Bool_t L1SingleMu18() const {return L1SingleMu18_;}
  Bool_t L2Mu10() const {return L2Mu10_;}
  Bool_t L1SingleMuOpen() const {return L1SingleMuOpen_;}
  Bool_t L1SingleMuOpen() const {return L1SingleMuOpen_;}
  Bool_t L2DoubleMu23() const {return L2DoubleMu23_;}
  Bool_t L2DoubleMu28() const {return L2DoubleMu28_;}
  Bool_t L2DoubleMu38() const {return L2DoubleMu38_;}
  Bool_t L2Mu10() const {return L2Mu10_;}
  Bool_t L2Mu10() const {return L2Mu10_;}
  Bool_t L2Mu45() const {return L2Mu45_;}
  Bool_t L2Mu40() const {return L2Mu40_;}
  Bool_t LooseIsoPFTau50() const {return LooseIsoPFTau50_;}
  Bool_t LooseIsoPFTau50() const {return LooseIsoPFTau50_;}
  Bool_t LooseIsoPFTau50() const {return LooseIsoPFTau50_;}
  Bool_t LooseIsoPFTau50() const {return LooseIsoPFTau50_;}
  Bool_t LooseIsoPFTau50() const {return LooseIsoPFTau50_;}
  Bool_t PFTau120() const {return PFTau120_;}
  Bool_t PFTau140() const {return PFTau140_;}
  Bool_t VLooseIsoPFTau120() const {return VLooseIsoPFTau120_;}
  Bool_t VLooseIsoPFTau140() const {return VLooseIsoPFTau140_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu20() const {return Mu20_;}
  Bool_t Mu20() const {return Mu20_;}
  Bool_t Mu20() const {return Mu20_;}
  Bool_t Mu20() const {return Mu20_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu25() const {return Mu25_;}
  Bool_t Mu27() const {return Mu27_;}
  Bool_t Mu30() const {return Mu30_;}
  Bool_t Mu30() const {return Mu30_;}
  Bool_t Mu40() const {return Mu40_;}
  Bool_t Mu40() const {return Mu40_;}
  Bool_t Mu20() const {return Mu20_;}
  Bool_t TkMu17() const {return TkMu17_;}
  Bool_t TkMu17() const {return TkMu17_;}
  Bool_t TkMu17() const {return TkMu17_;}
  Bool_t TkMu20() const {return TkMu20_;}
  Bool_t Mu24() const {return Mu24_;}
  Bool_t TkMu24() const {return TkMu24_;}
  Bool_t Mu27() const {return Mu27_;}
  Bool_t TkMu27() const {return TkMu27_;}
  Bool_t Mu45() const {return Mu45_;}
  Bool_t Mu50() const {return Mu50_;}
  Bool_t TkMu50() const {return TkMu50_;}
  Bool_t Mu38NoFiltersNoVtx() const {return Mu38NoFiltersNoVtx_;}
  Bool_t Mu42NoFiltersNoVtx() const {return Mu42NoFiltersNoVtx_;}
  Bool_t Mu28NoFiltersNoVtxDisplaced() const {return Mu28NoFiltersNoVtxDisplaced_;}
  Bool_t Mu33NoFiltersNoVtxDisplaced() const {return Mu33NoFiltersNoVtxDisplaced_;}
  Bool_t Mu23NoFiltersNoVtx() const {return Mu23NoFiltersNoVtx_;}
  Bool_t DoubleMu18NoFiltersNoVtx() const {return DoubleMu18NoFiltersNoVtx_;}
  Bool_t Mu33NoFiltersNoVtxDisplaced() const {return Mu33NoFiltersNoVtxDisplaced_;}
  Bool_t Mu33NoFiltersNoVtxDisplaced() const {return Mu33NoFiltersNoVtxDisplaced_;}
  Bool_t Mu28NoFiltersNoVtx() const {return Mu28NoFiltersNoVtx_;}
  Bool_t Mu38NoFiltersNoVtxDisplaced() const {return Mu38NoFiltersNoVtxDisplaced_;}
  Bool_t Mu38NoFiltersNoVtxDisplaced() const {return Mu38NoFiltersNoVtxDisplaced_;}
  Bool_t Mu38NoFiltersNoVtx() const {return Mu38NoFiltersNoVtx_;}
  Bool_t Mu28NoFiltersNoVtx() const {return Mu28NoFiltersNoVtx_;}
  Bool_t PFHT300() const {return PFHT300_;}
  Bool_t PFHT300() const {return PFHT300_;}
  Bool_t PFHT550() const {return PFHT550_;}
  Bool_t PFHT650() const {return PFHT650_;}
  Bool_t PFHT750() const {return PFHT750_;}
  Bool_t PFHT750() const {return PFHT750_;}
  Bool_t PFHT750() const {return PFHT750_;}
  Bool_t PFHT800() const {return PFHT800_;}
  Bool_t PFHT850() const {return PFHT850_;}
  Bool_t PFJet15() const {return PFJet15_;}
  Bool_t PFJet25() const {return PFJet25_;}
  Bool_t DiPFJet15() const {return DiPFJet15_;}
  Bool_t DiPFJet25() const {return DiPFJet25_;}
  Bool_t DiPFJet15() const {return DiPFJet15_;}
  Bool_t DiPFJet25() const {return DiPFJet25_;}
  Bool_t DiPFJetAve15() const {return DiPFJetAve15_;}
  Bool_t DiPFJetAve25() const {return DiPFJetAve25_;}
  Bool_t DiPFJetAve35() const {return DiPFJetAve35_;}
  Bool_t AK8PFJet40() const {return AK8PFJet40_;}
  Bool_t AK8PFJet60() const {return AK8PFJet60_;}
  Bool_t AK8PFJet80() const {return AK8PFJet80_;}
  Bool_t AK8PFJet140() const {return AK8PFJet140_;}
  Bool_t AK8PFJet200() const {return AK8PFJet200_;}
  Bool_t AK8PFJet260() const {return AK8PFJet260_;}
  Bool_t AK8PFJet320() const {return AK8PFJet320_;}
  Bool_t AK8PFJet400() const {return AK8PFJet400_;}
  Bool_t AK8PFJet450() const {return AK8PFJet450_;}
  Bool_t AK8PFJet500() const {return AK8PFJet500_;}
  Bool_t PFJet40() const {return PFJet40_;}
  Bool_t PFJet60() const {return PFJet60_;}
  Bool_t PFJet80() const {return PFJet80_;}
  Bool_t PFJet140() const {return PFJet140_;}
  Bool_t PFJet200() const {return PFJet200_;}
  Bool_t PFJet260() const {return PFJet260_;}
  Bool_t PFJet320() const {return PFJet320_;}
  Bool_t PFJet400() const {return PFJet400_;}
  Bool_t PFJet450() const {return PFJet450_;}
  Bool_t PFJet500() const {return PFJet500_;}
  Bool_t DiPFJetAve40() const {return DiPFJetAve40_;}
  Bool_t DiPFJetAve60() const {return DiPFJetAve60_;}
  Bool_t DiPFJetAve80() const {return DiPFJetAve80_;}
  Bool_t DiPFJetAve140() const {return DiPFJetAve140_;}
  Bool_t DiPFJetAve200() const {return DiPFJetAve200_;}
  Bool_t DiPFJetAve260() const {return DiPFJetAve260_;}
  Bool_t DiPFJetAve320() const {return DiPFJetAve320_;}
  Bool_t DiPFJetAve400() const {return DiPFJetAve400_;}
  Bool_t DiPFJetAve500() const {return DiPFJetAve500_;}
  Bool_t DiPFJetAve60() const {return DiPFJetAve60_;}
  Bool_t DiPFJetAve80() const {return DiPFJetAve80_;}
  Bool_t DiPFJetAve100() const {return DiPFJetAve100_;}
  Bool_t DiPFJetAve160() const {return DiPFJetAve160_;}
  Bool_t DiPFJetAve220() const {return DiPFJetAve220_;}
  Bool_t DiPFJetAve300() const {return DiPFJetAve300_;}
  Bool_t DiPFJet40() const {return DiPFJet40_;}
  Bool_t DiPFJet40() const {return DiPFJet40_;}
  Bool_t DiCentralPFJet170() const {return DiCentralPFJet170_;}
  Bool_t SingleCentralPFJet170() const {return SingleCentralPFJet170_;}
  Bool_t DiCentralPFJet170() const {return DiCentralPFJet170_;}
  Bool_t DiCentralPFJet220() const {return DiCentralPFJet220_;}
  Bool_t DiCentralPFJet330() const {return DiCentralPFJet330_;}
  Bool_t DiCentralPFJet430() const {return DiCentralPFJet430_;}
  Bool_t PFHT125() const {return PFHT125_;}
  Bool_t PFHT200() const {return PFHT200_;}
  Bool_t PFHT250() const {return PFHT250_;}
  Bool_t PFHT300() const {return PFHT300_;}
  Bool_t PFHT350() const {return PFHT350_;}
  Bool_t PFHT400() const {return PFHT400_;}
  Bool_t PFHT475() const {return PFHT475_;}
  Bool_t PFHT600() const {return PFHT600_;}
  Bool_t PFHT650() const {return PFHT650_;}
  Bool_t PFHT800() const {return PFHT800_;}
  Bool_t PFHT900() const {return PFHT900_;}
  Bool_t PFHT200() const {return PFHT200_;}
  Bool_t PFHT200() const {return PFHT200_;}
  Bool_t PFHT200() const {return PFHT200_;}
  Bool_t PFHT250() const {return PFHT250_;}
  Bool_t PFHT250() const {return PFHT250_;}
  Bool_t PFHT300() const {return PFHT300_;}
  Bool_t PFHT300() const {return PFHT300_;}
  Bool_t PFHT350() const {return PFHT350_;}
  Bool_t PFHT350() const {return PFHT350_;}
  Bool_t PFHT400() const {return PFHT400_;}
  Bool_t PFHT400() const {return PFHT400_;}
  Bool_t MET60() const {return MET60_;}
  Bool_t MET75() const {return MET75_;}
  Bool_t MET90() const {return MET90_;}
  Bool_t PFMET120() const {return PFMET120_;}
  Bool_t PFMET120() const {return PFMET120_;}
  Bool_t PFMET170() const {return PFMET170_;}
  Bool_t PFMET170() const {return PFMET170_;}
  Bool_t PFMET170() const {return PFMET170_;}
  Bool_t PFMET170() const {return PFMET170_;}
  Bool_t PFMET170() const {return PFMET170_;}
  Bool_t PFMET170() const {return PFMET170_;}
  Bool_t PFMETTypeOne190() const {return PFMETTypeOne190_;}
  Bool_t PFMET90() const {return PFMET90_;}
  Bool_t PFMET100() const {return PFMET100_;}
  Bool_t PFMET100() const {return PFMET100_;}
  Bool_t PFMET110() const {return PFMET110_;}
  Bool_t PFMET120() const {return PFMET120_;}
  Bool_t CaloMHTNoPU90() const {return CaloMHTNoPU90_;}
  Bool_t CaloMHTNoPU90() const {return CaloMHTNoPU90_;}
  Bool_t QuadPFJet() const {return QuadPFJet_;}
  Bool_t QuadPFJet() const {return QuadPFJet_;}
  Bool_t QuadPFJet() const {return QuadPFJet_;}
  Bool_t QuadPFJet() const {return QuadPFJet_;}
  Bool_t QuadPFJet() const {return QuadPFJet_;}
  Bool_t L1() const {return L1_;}
  Bool_t QuadJet45() const {return QuadJet45_;}
  Bool_t QuadJet45() const {return QuadJet45_;}
  Bool_t DoubleJet90() const {return DoubleJet90_;}
  Bool_t DoubleJet90() const {return DoubleJet90_;}
  Bool_t DoubleJetsC100() const {return DoubleJetsC100_;}
  Bool_t DoubleJetsC100() const {return DoubleJetsC100_;}
  Bool_t DoubleJetsC112() const {return DoubleJetsC112_;}
  Bool_t DoubleJetsC112() const {return DoubleJetsC112_;}
  Bool_t DoubleJetsC100() const {return DoubleJetsC100_;}
  Bool_t DoubleJetsC100() const {return DoubleJetsC100_;}
  Bool_t DoubleJetsC100() const {return DoubleJetsC100_;}
  Bool_t DoubleJetsC100() const {return DoubleJetsC100_;}
  Bool_t Photon135() const {return Photon135_;}
  Bool_t Photon20() const {return Photon20_;}
  Bool_t Photon22() const {return Photon22_;}
  Bool_t Photon22() const {return Photon22_;}
  Bool_t Photon250() const {return Photon250_;}
  Bool_t Photon300() const {return Photon300_;}
  Bool_t Photon26() const {return Photon26_;}
  Bool_t Photon36() const {return Photon36_;}
  Bool_t Photon36() const {return Photon36_;}
  Bool_t Photon36() const {return Photon36_;}
  Bool_t Photon50() const {return Photon50_;}
  Bool_t Photon50() const {return Photon50_;}
  Bool_t Photon75() const {return Photon75_;}
  Bool_t Photon75() const {return Photon75_;}
  Bool_t Photon90() const {return Photon90_;}
  Bool_t Photon90() const {return Photon90_;}
  Bool_t Photon120() const {return Photon120_;}
  Bool_t Photon120() const {return Photon120_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Ele8() const {return Ele8_;}
  Bool_t Ele12() const {return Ele12_;}
  Bool_t Ele17() const {return Ele17_;}
  Bool_t Ele23() const {return Ele23_;}
  Bool_t BTagMu() const {return BTagMu_;}
  Bool_t BTagMu() const {return BTagMu_;}
  Bool_t BTagMu() const {return BTagMu_;}
  Bool_t BTagMu() const {return BTagMu_;}
  Bool_t BTagMu() const {return BTagMu_;}
  Bool_t BTagMu() const {return BTagMu_;}
  Bool_t BTagMu() const {return BTagMu_;}
  Bool_t Ele23() const {return Ele23_;}
  Bool_t Ele23() const {return Ele23_;}
  Bool_t Ele17() const {return Ele17_;}
  Bool_t Ele16() const {return Ele16_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t Mu12() const {return Mu12_;}
  Bool_t Mu12() const {return Mu12_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu23() const {return Mu23_;}
  Bool_t Mu23() const {return Mu23_;}
  Bool_t Mu23() const {return Mu23_;}
  Bool_t Mu23() const {return Mu23_;}
  Bool_t Mu30() const {return Mu30_;}
  Bool_t Mu33() const {return Mu33_;}
  Bool_t Mu37() const {return Mu37_;}
  Bool_t Mu27() const {return Mu27_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t Mu12() const {return Mu12_;}
  Bool_t Mu12() const {return Mu12_;}
  Bool_t Mu12() const {return Mu12_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t DiMu9() const {return DiMu9_;}
  Bool_t TripleMu() const {return TripleMu_;}
  Bool_t TripleMu() const {return TripleMu_;}
  Bool_t Mu3er() const {return Mu3er_;}
  Bool_t Mu6() const {return Mu6_;}
  Bool_t Mu6() const {return Mu6_;}
  Bool_t Mu14er() const {return Mu14er_;}
  Bool_t Ele17() const {return Ele17_;}
  Bool_t Ele23() const {return Ele23_;}
  Bool_t Ele12() const {return Ele12_;}
  Bool_t Ele17() const {return Ele17_;}
  Bool_t Ele17() const {return Ele17_;}
  Bool_t Ele23() const {return Ele23_;}
  Bool_t PFHT650() const {return PFHT650_;}
  Bool_t PFHT650() const {return PFHT650_;}
  Bool_t Photon22() const {return Photon22_;}
  Bool_t Photon30() const {return Photon30_;}
  Bool_t Photon36() const {return Photon36_;}
  Bool_t Photon50() const {return Photon50_;}
  Bool_t Photon75() const {return Photon75_;}
  Bool_t Photon90() const {return Photon90_;}
  Bool_t Photon120() const {return Photon120_;}
  Bool_t Photon175() const {return Photon175_;}
  Bool_t Photon165() const {return Photon165_;}
  Bool_t Photon22() const {return Photon22_;}
  Bool_t Photon30() const {return Photon30_;}
  Bool_t Photon36() const {return Photon36_;}
  Bool_t Photon50() const {return Photon50_;}
  Bool_t Photon75() const {return Photon75_;}
  Bool_t Photon90() const {return Photon90_;}
  Bool_t Photon120() const {return Photon120_;}
  Bool_t Photon165() const {return Photon165_;}
  Bool_t Diphoton30() const {return Diphoton30_;}
  Bool_t Diphoton30() const {return Diphoton30_;}
  Bool_t Diphoton30PV() const {return Diphoton30PV_;}
  Bool_t Diphoton30() const {return Diphoton30_;}
  Bool_t Diphoton30EB() const {return Diphoton30EB_;}
  Bool_t Dimuon0() const {return Dimuon0_;}
  Bool_t Dimuon0() const {return Dimuon0_;}
  Bool_t QuadMuon0() const {return QuadMuon0_;}
  Bool_t QuadMuon0() const {return QuadMuon0_;}
  Bool_t Rsq0p25() const {return Rsq0p25_;}
  Bool_t RsqMR240() const {return RsqMR240_;}
  Bool_t RsqMR240() const {return RsqMR240_;}
  Bool_t Rsq0p25() const {return Rsq0p25_;}
  Bool_t Rsq0p30() const {return Rsq0p30_;}
  Bool_t RsqMR240() const {return RsqMR240_;}
  Bool_t RsqMR240() const {return RsqMR240_;}
  Bool_t RsqMR270() const {return RsqMR270_;}
  Bool_t RsqMR270() const {return RsqMR270_;}
  Bool_t Rsq0p02() const {return Rsq0p02_;}
  Bool_t Rsq0p02() const {return Rsq0p02_;}
  Bool_t Rsq0p02() const {return Rsq0p02_;}
  Bool_t Rsq0p02() const {return Rsq0p02_;}
  Bool_t Rsq0p02() const {return Rsq0p02_;}
  Bool_t HT200() const {return HT200_;}
  Bool_t HT250() const {return HT250_;}
  Bool_t HT350() const {return HT350_;}
  Bool_t HT350() const {return HT350_;}
  Bool_t HT350() const {return HT350_;}
  Bool_t HT350() const {return HT350_;}
  Bool_t HT400() const {return HT400_;}
  Bool_t HT500() const {return HT500_;}
  Bool_t HT550() const {return HT550_;}
  Bool_t HT550() const {return HT550_;}
  Bool_t HT650() const {return HT650_;}
  Bool_t HT750() const {return HT750_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t VBF() const {return VBF_;}
  Bool_t PFMETNoMu90() const {return PFMETNoMu90_;}
  Bool_t PFMETNoMu100() const {return PFMETNoMu100_;}
  Bool_t PFMETNoMu110() const {return PFMETNoMu110_;}
  Bool_t PFMETNoMu120() const {return PFMETNoMu120_;}
  Bool_t MonoCentralPFJet80() const {return MonoCentralPFJet80_;}
  Bool_t MonoCentralPFJet80() const {return MonoCentralPFJet80_;}
  Bool_t MonoCentralPFJet80() const {return MonoCentralPFJet80_;}
  Bool_t MonoCentralPFJet80() const {return MonoCentralPFJet80_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t Photon90() const {return Photon90_;}
  Bool_t DoubleMu8() const {return DoubleMu8_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t DoubleEle8() const {return DoubleEle8_;}
  Bool_t DoubleMu8() const {return DoubleMu8_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t DoubleEle8() const {return DoubleEle8_;}
  Bool_t Mu10() const {return Mu10_;}
  Bool_t DoubleMu3() const {return DoubleMu3_;}
  Bool_t Ele10() const {return Ele10_;}
  Bool_t Ele15() const {return Ele15_;}
  Bool_t Ele15() const {return Ele15_;}
  Bool_t Ele15() const {return Ele15_;}
  Bool_t Ele15() const {return Ele15_;}
  Bool_t Ele15() const {return Ele15_;}
  Bool_t Ele15() const {return Ele15_;}
  Bool_t Ele50() const {return Ele50_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t Mu10() const {return Mu10_;}
  Bool_t Mu15() const {return Mu15_;}
  Bool_t Mu15() const {return Mu15_;}
  Bool_t Mu15() const {return Mu15_;}
  Bool_t Mu15() const {return Mu15_;}
  Bool_t Mu15() const {return Mu15_;}
  Bool_t Mu15() const {return Mu15_;}
  Bool_t Mu50() const {return Mu50_;}
  Bool_t Dimuon16() const {return Dimuon16_;}
  Bool_t Dimuon10() const {return Dimuon10_;}
  Bool_t Dimuon8() const {return Dimuon8_;}
  Bool_t Dimuon8() const {return Dimuon8_;}
  Bool_t Dimuon0() const {return Dimuon0_;}
  Bool_t Mu16() const {return Mu16_;}
  Bool_t Mu16() const {return Mu16_;}
  Bool_t TrkMu15() const {return TrkMu15_;}
  Bool_t TrkMu17() const {return TrkMu17_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu3() const {return Mu3_;}
  Bool_t Ele8() const {return Ele8_;}
  Bool_t Ele12() const {return Ele12_;}
  Bool_t Ele17() const {return Ele17_;}
  Bool_t Ele23() const {return Ele23_;}
  Bool_t Ele50() const {return Ele50_;}
  Bool_t Ele50() const {return Ele50_;}
  Bool_t PFHT400() const {return PFHT400_;}
  Bool_t PFHT450() const {return PFHT450_;}
  Bool_t PFHT400() const {return PFHT400_;}
  Bool_t PFHT450() const {return PFHT450_;}
  Bool_t Ele115() const {return Ele115_;}
  Bool_t Mu55() const {return Mu55_;}
  Bool_t Photon42() const {return Photon42_;}
  Bool_t Photon90() const {return Photon90_;}
  Bool_t PixelTracks() const {return PixelTracks_;}
  Bool_t PixelTracks() const {return PixelTracks_;}
  Bool_t PixelTracks() const {return PixelTracks_;}
  Bool_t PixelTracks() const {return PixelTracks_;}
  Bool_t PixelTracks() const {return PixelTracks_;}
  Bool_t FullTracks() const {return FullTracks_;}
  Bool_t FullTracks() const {return FullTracks_;}
  Bool_t FullTracks() const {return FullTracks_;}
  Bool_t FullTracks() const {return FullTracks_;}
  Bool_t ECALHT800() const {return ECALHT800_;}
  Bool_t DiSC30() const {return DiSC30_;}
  Bool_t Photon125() const {return Photon125_;}
  Bool_t MET100() const {return MET100_;}
  Bool_t MET150() const {return MET150_;}
  Bool_t MET200() const {return MET200_;}
  Bool_t Ele27() const {return Ele27_;}
  Bool_t L1FatEvents() const {return L1FatEvents_;}
  Bool_t Physics() const {return Physics_;}
  Bool_t L1FatEvents() const {return L1FatEvents_;}
  Bool_t L1FatEvents() const {return L1FatEvents_;}
  Bool_t L1FatEvents() const {return L1FatEvents_;}
  Bool_t L1FatEvents() const {return L1FatEvents_;}
  Bool_t Random() const {return Random_;}
  Bool_t ZeroBias() const {return ZeroBias_;}
  Bool_t AK4CaloJet30() const {return AK4CaloJet30_;}
  Bool_t AK4CaloJet40() const {return AK4CaloJet40_;}
  Bool_t AK4CaloJet50() const {return AK4CaloJet50_;}
  Bool_t AK4CaloJet80() const {return AK4CaloJet80_;}
  Bool_t AK4CaloJet100() const {return AK4CaloJet100_;}
  Bool_t AK4PFJet30() const {return AK4PFJet30_;}
  Bool_t AK4PFJet50() const {return AK4PFJet50_;}
  Bool_t AK4PFJet80() const {return AK4PFJet80_;}
  Bool_t AK4PFJet100() const {return AK4PFJet100_;}
  Bool_t HISinglePhoton10() const {return HISinglePhoton10_;}
  Bool_t HISinglePhoton15() const {return HISinglePhoton15_;}
  Bool_t HISinglePhoton20() const {return HISinglePhoton20_;}
  Bool_t HISinglePhoton40() const {return HISinglePhoton40_;}
  Bool_t HISinglePhoton60() const {return HISinglePhoton60_;}
  Bool_t EcalCalibration() const {return EcalCalibration_;}
  Bool_t HcalCalibration() const {return HcalCalibration_;}
  Bool_t GlobalRunHPDNoise() const {return GlobalRunHPDNoise_;}
  Bool_t L1BptxMinus() const {return L1BptxMinus_;}
  Bool_t L1BptxPlus() const {return L1BptxPlus_;}
  Bool_t L1NotBptxOR() const {return L1NotBptxOR_;}
  Bool_t L1BeamGasMinus() const {return L1BeamGasMinus_;}
  Bool_t L1BeamGasPlus() const {return L1BeamGasPlus_;}
  Bool_t L1BptxXOR() const {return L1BptxXOR_;}
  Bool_t L1MinimumBiasHF() const {return L1MinimumBiasHF_;}
  Bool_t L1MinimumBiasHF() const {return L1MinimumBiasHF_;}
  Bool_t HcalNZS() const {return HcalNZS_;}
  Bool_t HcalPhiSym() const {return HcalPhiSym_;}
  Bool_t HcalIsolatedbunch() const {return HcalIsolatedbunch_;}
  Bool_t ZeroBias() const {return ZeroBias_;}
  Bool_t ZeroBias() const {return ZeroBias_;}
  Bool_t ZeroBias() const {return ZeroBias_;}
  Bool_t ZeroBias() const {return ZeroBias_;}
  Bool_t ZeroBias() const {return ZeroBias_;}
  Bool_t ZeroBias() const {return ZeroBias_;}
  Bool_t Photon500() const {return Photon500_;}
  Bool_t Photon600() const {return Photon600_;}
  Bool_t Mu300() const {return Mu300_;}
  Bool_t Mu350() const {return Mu350_;}
  Bool_t MET250() const {return MET250_;}
  Bool_t MET300() const {return MET300_;}
  Bool_t MET600() const {return MET600_;}
  Bool_t MET700() const {return MET700_;}
  Bool_t PFMET300() const {return PFMET300_;}
  Bool_t PFMET400() const {return PFMET400_;}
  Bool_t PFMET500() const {return PFMET500_;}
  Bool_t PFMET600() const {return PFMET600_;}
  Bool_t Ele250() const {return Ele250_;}
  Bool_t Ele300() const {return Ele300_;}
  Bool_t HT2000() const {return HT2000_;}
  Bool_t HT2500() const {return HT2500_;}
  Bool_t IsoTrackHE() const {return IsoTrackHE_;}
  Bool_t IsoTrackHB() const {return IsoTrackHB_;}
  Bool_t HLTriggerFinalPath() const {return HLTriggerFinalPath_;}
private:
  Bool_t HLTriggerFirstPath_;
  Bool_t AK8PFJet360_;
  Bool_t AK8PFJet400_;
  Bool_t AK8PFHT750_;
  Bool_t AK8PFHT800_;
  Bool_t AK8DiPFJet300_;
  Bool_t AK8DiPFJet280_;
  Bool_t AK8DiPFJet300_;
  Bool_t AK8DiPFJet300_;
  Bool_t AK8PFHT700_;
  Bool_t AK8PFHT650_;
  Bool_t AK8PFHT600_;
  Bool_t AK8DiPFJet280_;
  Bool_t AK8DiPFJet250_;
  Bool_t AK8DiPFJet280_;
  Bool_t AK8DiPFJet250_;
  Bool_t CaloJet260_;
  Bool_t CaloJet500_;
  Bool_t Dimuon13_;
  Bool_t Dimuon13_;
  Bool_t Dimuon20_;
  Bool_t DoubleEle24_;
  Bool_t DoubleEle25_;
  Bool_t DoubleEle33_;
  Bool_t DoubleEle33_;
  Bool_t DoubleEle33_;
  Bool_t DoubleEle33_;
  Bool_t DoubleMediumCombinedIsoPFTau35_;
  Bool_t DoubleTightCombinedIsoPFTau35_;
  Bool_t DoubleMediumCombinedIsoPFTau40_;
  Bool_t DoubleTightCombinedIsoPFTau40_;
  Bool_t DoubleMediumCombinedIsoPFTau40_;
  Bool_t DoubleTightCombinedIsoPFTau40_;
  Bool_t DoubleMediumIsoPFTau35_;
  Bool_t DoubleMediumIsoPFTau40_;
  Bool_t DoubleMediumIsoPFTau40_;
  Bool_t DoubleEle37_;
  Bool_t DoubleMu33NoFiltersNoVtx_;
  Bool_t DoubleMu38NoFiltersNoVtx_;
  Bool_t DoubleMu23NoFiltersNoVtxDisplaced_;
  Bool_t DoubleMu28NoFiltersNoVtxDisplaced_;
  Bool_t DoubleMu0_;
  Bool_t DoubleMu4_;
  Bool_t DoubleMu4_;
  Bool_t DoubleMu4_;
  Bool_t DoubleMu4_;
  Bool_t DoubleMu3_;
  Bool_t DoubleMu4_;
  Bool_t Mu7p5_;
  Bool_t Mu7p5_;
  Bool_t Mu7p5_;
  Bool_t Mu7p5_;
  Bool_t Mu7p5_;
  Bool_t Mu7p5_;
  Bool_t Mu7p5_;
  Bool_t Mu7p5_;
  Bool_t Dimuon0er16_;
  Bool_t Dimuon0er16_;
  Bool_t Dimuon6_;
  Bool_t Photon150_;
  Bool_t Photon90_;
  Bool_t HT250_;
  Bool_t DoublePhoton60_;
  Bool_t DoublePhoton85_;
  Bool_t Ele17_;
  Bool_t Ele20_;
  Bool_t Ele22_;
  Bool_t Ele22_;
  Bool_t Ele22_;
  Bool_t Ele23_;
  Bool_t Ele23_;
  Bool_t Ele24_;
  Bool_t Ele24_;
  Bool_t Ele24_;
  Bool_t Ele24_;
  Bool_t Ele25_;
  Bool_t Ele25_;
  Bool_t Ele25_;
  Bool_t Ele27_;
  Bool_t Ele27_;
  Bool_t Ele27_;
  Bool_t Ele27_;
  Bool_t Ele27_;
  Bool_t Ele27_;
  Bool_t Ele27_;
  Bool_t Ele30_;
  Bool_t Ele30_;
  Bool_t Ele30_;
  Bool_t Ele32_;
  Bool_t Ele32_;
  Bool_t Ele32_;
  Bool_t Ele32_;
  Bool_t Ele35_;
  Bool_t Ele35_;
  Bool_t Ele36_;
  Bool_t Ele45_;
  Bool_t Ele45_;
  Bool_t Ele45_;
  Bool_t Ele105_;
  Bool_t Ele30WP60_;
  Bool_t Ele30WP60_;
  Bool_t HT200_;
  Bool_t HT275_;
  Bool_t HT325_;
  Bool_t HT425_;
  Bool_t HT575_;
  Bool_t HT410to430_;
  Bool_t HT430to450_;
  Bool_t HT450to470_;
  Bool_t HT470to500_;
  Bool_t HT500to550_;
  Bool_t HT550to650_;
  Bool_t HT650_;
  Bool_t Mu16_;
  Bool_t IsoMu16_;
  Bool_t IsoMu16_;
  Bool_t IsoMu17_;
  Bool_t IsoMu17_;
  Bool_t IsoMu17_;
  Bool_t DoubleIsoMu17_;
  Bool_t DoubleIsoMu17_;
  Bool_t IsoMu18_;
  Bool_t IsoMu19_;
  Bool_t IsoMu19_;
  Bool_t IsoMu19_;
  Bool_t IsoMu19_;
  Bool_t IsoMu19_;
  Bool_t IsoMu19_;
  Bool_t IsoMu21_;
  Bool_t IsoMu21_;
  Bool_t IsoMu20_;
  Bool_t IsoMu21_;
  Bool_t IsoMu21_;
  Bool_t IsoMu21_;
  Bool_t IsoMu22_;
  Bool_t IsoMu22_;
  Bool_t IsoMu24_;
  Bool_t IsoMu27_;
  Bool_t IsoTkMu18_;
  Bool_t IsoTkMu20_;
  Bool_t IsoTkMu22_;
  Bool_t IsoTkMu22_;
  Bool_t IsoTkMu24_;
  Bool_t IsoTkMu27_;
  Bool_t JetE30_;
  Bool_t JetE30_;
  Bool_t JetE50_;
  Bool_t JetE70_;
  Bool_t L1SingleMu18_;
  Bool_t L2Mu10_;
  Bool_t L1SingleMuOpen_;
  Bool_t L1SingleMuOpen_;
  Bool_t L2DoubleMu23_;
  Bool_t L2DoubleMu28_;
  Bool_t L2DoubleMu38_;
  Bool_t L2Mu10_;
  Bool_t L2Mu10_;
  Bool_t L2Mu45_;
  Bool_t L2Mu40_;
  Bool_t LooseIsoPFTau50_;
  Bool_t LooseIsoPFTau50_;
  Bool_t LooseIsoPFTau50_;
  Bool_t LooseIsoPFTau50_;
  Bool_t LooseIsoPFTau50_;
  Bool_t PFTau120_;
  Bool_t PFTau140_;
  Bool_t VLooseIsoPFTau120_;
  Bool_t VLooseIsoPFTau140_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t Mu20_;
  Bool_t Mu20_;
  Bool_t Mu20_;
  Bool_t Mu20_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t Mu25_;
  Bool_t Mu27_;
  Bool_t Mu30_;
  Bool_t Mu30_;
  Bool_t Mu40_;
  Bool_t Mu40_;
  Bool_t Mu20_;
  Bool_t TkMu17_;
  Bool_t TkMu17_;
  Bool_t TkMu17_;
  Bool_t TkMu20_;
  Bool_t Mu24_;
  Bool_t TkMu24_;
  Bool_t Mu27_;
  Bool_t TkMu27_;
  Bool_t Mu45_;
  Bool_t Mu50_;
  Bool_t TkMu50_;
  Bool_t Mu38NoFiltersNoVtx_;
  Bool_t Mu42NoFiltersNoVtx_;
  Bool_t Mu28NoFiltersNoVtxDisplaced_;
  Bool_t Mu33NoFiltersNoVtxDisplaced_;
  Bool_t Mu23NoFiltersNoVtx_;
  Bool_t DoubleMu18NoFiltersNoVtx_;
  Bool_t Mu33NoFiltersNoVtxDisplaced_;
  Bool_t Mu33NoFiltersNoVtxDisplaced_;
  Bool_t Mu28NoFiltersNoVtx_;
  Bool_t Mu38NoFiltersNoVtxDisplaced_;
  Bool_t Mu38NoFiltersNoVtxDisplaced_;
  Bool_t Mu38NoFiltersNoVtx_;
  Bool_t Mu28NoFiltersNoVtx_;
  Bool_t PFHT300_;
  Bool_t PFHT300_;
  Bool_t PFHT550_;
  Bool_t PFHT650_;
  Bool_t PFHT750_;
  Bool_t PFHT750_;
  Bool_t PFHT750_;
  Bool_t PFHT800_;
  Bool_t PFHT850_;
  Bool_t PFJet15_;
  Bool_t PFJet25_;
  Bool_t DiPFJet15_;
  Bool_t DiPFJet25_;
  Bool_t DiPFJet15_;
  Bool_t DiPFJet25_;
  Bool_t DiPFJetAve15_;
  Bool_t DiPFJetAve25_;
  Bool_t DiPFJetAve35_;
  Bool_t AK8PFJet40_;
  Bool_t AK8PFJet60_;
  Bool_t AK8PFJet80_;
  Bool_t AK8PFJet140_;
  Bool_t AK8PFJet200_;
  Bool_t AK8PFJet260_;
  Bool_t AK8PFJet320_;
  Bool_t AK8PFJet400_;
  Bool_t AK8PFJet450_;
  Bool_t AK8PFJet500_;
  Bool_t PFJet40_;
  Bool_t PFJet60_;
  Bool_t PFJet80_;
  Bool_t PFJet140_;
  Bool_t PFJet200_;
  Bool_t PFJet260_;
  Bool_t PFJet320_;
  Bool_t PFJet400_;
  Bool_t PFJet450_;
  Bool_t PFJet500_;
  Bool_t DiPFJetAve40_;
  Bool_t DiPFJetAve60_;
  Bool_t DiPFJetAve80_;
  Bool_t DiPFJetAve140_;
  Bool_t DiPFJetAve200_;
  Bool_t DiPFJetAve260_;
  Bool_t DiPFJetAve320_;
  Bool_t DiPFJetAve400_;
  Bool_t DiPFJetAve500_;
  Bool_t DiPFJetAve60_;
  Bool_t DiPFJetAve80_;
  Bool_t DiPFJetAve100_;
  Bool_t DiPFJetAve160_;
  Bool_t DiPFJetAve220_;
  Bool_t DiPFJetAve300_;
  Bool_t DiPFJet40_;
  Bool_t DiPFJet40_;
  Bool_t DiCentralPFJet170_;
  Bool_t SingleCentralPFJet170_;
  Bool_t DiCentralPFJet170_;
  Bool_t DiCentralPFJet220_;
  Bool_t DiCentralPFJet330_;
  Bool_t DiCentralPFJet430_;
  Bool_t PFHT125_;
  Bool_t PFHT200_;
  Bool_t PFHT250_;
  Bool_t PFHT300_;
  Bool_t PFHT350_;
  Bool_t PFHT400_;
  Bool_t PFHT475_;
  Bool_t PFHT600_;
  Bool_t PFHT650_;
  Bool_t PFHT800_;
  Bool_t PFHT900_;
  Bool_t PFHT200_;
  Bool_t PFHT200_;
  Bool_t PFHT200_;
  Bool_t PFHT250_;
  Bool_t PFHT250_;
  Bool_t PFHT300_;
  Bool_t PFHT300_;
  Bool_t PFHT350_;
  Bool_t PFHT350_;
  Bool_t PFHT400_;
  Bool_t PFHT400_;
  Bool_t MET60_;
  Bool_t MET75_;
  Bool_t MET90_;
  Bool_t PFMET120_;
  Bool_t PFMET120_;
  Bool_t PFMET170_;
  Bool_t PFMET170_;
  Bool_t PFMET170_;
  Bool_t PFMET170_;
  Bool_t PFMET170_;
  Bool_t PFMET170_;
  Bool_t PFMETTypeOne190_;
  Bool_t PFMET90_;
  Bool_t PFMET100_;
  Bool_t PFMET100_;
  Bool_t PFMET110_;
  Bool_t PFMET120_;
  Bool_t CaloMHTNoPU90_;
  Bool_t CaloMHTNoPU90_;
  Bool_t QuadPFJet_;
  Bool_t QuadPFJet_;
  Bool_t QuadPFJet_;
  Bool_t QuadPFJet_;
  Bool_t QuadPFJet_;
  Bool_t L1_;
  Bool_t QuadJet45_;
  Bool_t QuadJet45_;
  Bool_t DoubleJet90_;
  Bool_t DoubleJet90_;
  Bool_t DoubleJetsC100_;
  Bool_t DoubleJetsC100_;
  Bool_t DoubleJetsC112_;
  Bool_t DoubleJetsC112_;
  Bool_t DoubleJetsC100_;
  Bool_t DoubleJetsC100_;
  Bool_t DoubleJetsC100_;
  Bool_t DoubleJetsC100_;
  Bool_t Photon135_;
  Bool_t Photon20_;
  Bool_t Photon22_;
  Bool_t Photon22_;
  Bool_t Photon250_;
  Bool_t Photon300_;
  Bool_t Photon26_;
  Bool_t Photon36_;
  Bool_t Photon36_;
  Bool_t Photon36_;
  Bool_t Photon50_;
  Bool_t Photon50_;
  Bool_t Photon75_;
  Bool_t Photon75_;
  Bool_t Photon90_;
  Bool_t Photon90_;
  Bool_t Photon120_;
  Bool_t Photon120_;
  Bool_t Mu8_;
  Bool_t Mu17_;
  Bool_t Ele8_;
  Bool_t Ele12_;
  Bool_t Ele17_;
  Bool_t Ele23_;
  Bool_t BTagMu_;
  Bool_t BTagMu_;
  Bool_t BTagMu_;
  Bool_t BTagMu_;
  Bool_t BTagMu_;
  Bool_t BTagMu_;
  Bool_t BTagMu_;
  Bool_t Ele23_;
  Bool_t Ele23_;
  Bool_t Ele17_;
  Bool_t Ele16_;
  Bool_t Mu8_;
  Bool_t Mu8_;
  Bool_t Mu8_;
  Bool_t Mu12_;
  Bool_t Mu12_;
  Bool_t Mu17_;
  Bool_t Mu23_;
  Bool_t Mu23_;
  Bool_t Mu23_;
  Bool_t Mu23_;
  Bool_t Mu30_;
  Bool_t Mu33_;
  Bool_t Mu37_;
  Bool_t Mu27_;
  Bool_t Mu8_;
  Bool_t Mu12_;
  Bool_t Mu12_;
  Bool_t Mu12_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t Mu17_;
  Bool_t DiMu9_;
  Bool_t TripleMu_;
  Bool_t TripleMu_;
  Bool_t Mu3er_;
  Bool_t Mu6_;
  Bool_t Mu6_;
  Bool_t Mu14er_;
  Bool_t Ele17_;
  Bool_t Ele23_;
  Bool_t Ele12_;
  Bool_t Ele17_;
  Bool_t Ele17_;
  Bool_t Ele23_;
  Bool_t PFHT650_;
  Bool_t PFHT650_;
  Bool_t Photon22_;
  Bool_t Photon30_;
  Bool_t Photon36_;
  Bool_t Photon50_;
  Bool_t Photon75_;
  Bool_t Photon90_;
  Bool_t Photon120_;
  Bool_t Photon175_;
  Bool_t Photon165_;
  Bool_t Photon22_;
  Bool_t Photon30_;
  Bool_t Photon36_;
  Bool_t Photon50_;
  Bool_t Photon75_;
  Bool_t Photon90_;
  Bool_t Photon120_;
  Bool_t Photon165_;
  Bool_t Diphoton30_;
  Bool_t Diphoton30_;
  Bool_t Diphoton30PV_;
  Bool_t Diphoton30_;
  Bool_t Diphoton30EB_;
  Bool_t Dimuon0_;
  Bool_t Dimuon0_;
  Bool_t QuadMuon0_;
  Bool_t QuadMuon0_;
  Bool_t Rsq0p25_;
  Bool_t RsqMR240_;
  Bool_t RsqMR240_;
  Bool_t Rsq0p25_;
  Bool_t Rsq0p30_;
  Bool_t RsqMR240_;
  Bool_t RsqMR240_;
  Bool_t RsqMR270_;
  Bool_t RsqMR270_;
  Bool_t Rsq0p02_;
  Bool_t Rsq0p02_;
  Bool_t Rsq0p02_;
  Bool_t Rsq0p02_;
  Bool_t Rsq0p02_;
  Bool_t HT200_;
  Bool_t HT250_;
  Bool_t HT350_;
  Bool_t HT350_;
  Bool_t HT350_;
  Bool_t HT350_;
  Bool_t HT400_;
  Bool_t HT500_;
  Bool_t HT550_;
  Bool_t HT550_;
  Bool_t HT650_;
  Bool_t HT750_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t VBF_;
  Bool_t PFMETNoMu90_;
  Bool_t PFMETNoMu100_;
  Bool_t PFMETNoMu110_;
  Bool_t PFMETNoMu120_;
  Bool_t MonoCentralPFJet80_;
  Bool_t MonoCentralPFJet80_;
  Bool_t MonoCentralPFJet80_;
  Bool_t MonoCentralPFJet80_;
  Bool_t Ele27_;
  Bool_t Photon90_;
  Bool_t DoubleMu8_;
  Bool_t Mu8_;
  Bool_t DoubleEle8_;
  Bool_t DoubleMu8_;
  Bool_t Mu8_;
  Bool_t DoubleEle8_;
  Bool_t Mu10_;
  Bool_t DoubleMu3_;
  Bool_t Ele10_;
  Bool_t Ele15_;
  Bool_t Ele15_;
  Bool_t Ele15_;
  Bool_t Ele15_;
  Bool_t Ele15_;
  Bool_t Ele15_;
  Bool_t Ele50_;
  Bool_t Mu8_;
  Bool_t Mu10_;
  Bool_t Mu15_;
  Bool_t Mu15_;
  Bool_t Mu15_;
  Bool_t Mu15_;
  Bool_t Mu15_;
  Bool_t Mu15_;
  Bool_t Mu50_;
  Bool_t Dimuon16_;
  Bool_t Dimuon10_;
  Bool_t Dimuon8_;
  Bool_t Dimuon8_;
  Bool_t Dimuon0_;
  Bool_t Mu16_;
  Bool_t Mu16_;
  Bool_t TrkMu15_;
  Bool_t TrkMu17_;
  Bool_t Mu8_;
  Bool_t Mu17_;
  Bool_t Mu3_;
  Bool_t Ele8_;
  Bool_t Ele12_;
  Bool_t Ele17_;
  Bool_t Ele23_;
  Bool_t Ele50_;
  Bool_t Ele50_;
  Bool_t PFHT400_;
  Bool_t PFHT450_;
  Bool_t PFHT400_;
  Bool_t PFHT450_;
  Bool_t Ele115_;
  Bool_t Mu55_;
  Bool_t Photon42_;
  Bool_t Photon90_;
  Bool_t PixelTracks_;
  Bool_t PixelTracks_;
  Bool_t PixelTracks_;
  Bool_t PixelTracks_;
  Bool_t PixelTracks_;
  Bool_t FullTracks_;
  Bool_t FullTracks_;
  Bool_t FullTracks_;
  Bool_t FullTracks_;
  Bool_t ECALHT800_;
  Bool_t DiSC30_;
  Bool_t Photon125_;
  Bool_t MET100_;
  Bool_t MET150_;
  Bool_t MET200_;
  Bool_t Ele27_;
  Bool_t L1FatEvents_;
  Bool_t Physics_;
  Bool_t L1FatEvents_;
  Bool_t L1FatEvents_;
  Bool_t L1FatEvents_;
  Bool_t L1FatEvents_;
  Bool_t Random_;
  Bool_t ZeroBias_;
  Bool_t AK4CaloJet30_;
  Bool_t AK4CaloJet40_;
  Bool_t AK4CaloJet50_;
  Bool_t AK4CaloJet80_;
  Bool_t AK4CaloJet100_;
  Bool_t AK4PFJet30_;
  Bool_t AK4PFJet50_;
  Bool_t AK4PFJet80_;
  Bool_t AK4PFJet100_;
  Bool_t HISinglePhoton10_;
  Bool_t HISinglePhoton15_;
  Bool_t HISinglePhoton20_;
  Bool_t HISinglePhoton40_;
  Bool_t HISinglePhoton60_;
  Bool_t EcalCalibration_;
  Bool_t HcalCalibration_;
  Bool_t GlobalRunHPDNoise_;
  Bool_t L1BptxMinus_;
  Bool_t L1BptxPlus_;
  Bool_t L1NotBptxOR_;
  Bool_t L1BeamGasMinus_;
  Bool_t L1BeamGasPlus_;
  Bool_t L1BptxXOR_;
  Bool_t L1MinimumBiasHF_;
  Bool_t L1MinimumBiasHF_;
  Bool_t HcalNZS_;
  Bool_t HcalPhiSym_;
  Bool_t HcalIsolatedbunch_;
  Bool_t ZeroBias_;
  Bool_t ZeroBias_;
  Bool_t ZeroBias_;
  Bool_t ZeroBias_;
  Bool_t ZeroBias_;
  Bool_t ZeroBias_;
  Bool_t Photon500_;
  Bool_t Photon600_;
  Bool_t Mu300_;
  Bool_t Mu350_;
  Bool_t MET250_;
  Bool_t MET300_;
  Bool_t MET600_;
  Bool_t MET700_;
  Bool_t PFMET300_;
  Bool_t PFMET400_;
  Bool_t PFMET500_;
  Bool_t PFMET600_;
  Bool_t Ele250_;
  Bool_t Ele300_;
  Bool_t HT2000_;
  Bool_t HT2500_;
  Bool_t IsoTrackHE_;
  Bool_t IsoTrackHB_;
  Bool_t HLTriggerFinalPath_;
  void setHLTriggerFirstPath(const Bool_t value) {HLTriggerFirstPath_ = value;}
  void setAK8PFJet360(const Bool_t value) {AK8PFJet360_ = value;}
  void setAK8PFJet400(const Bool_t value) {AK8PFJet400_ = value;}
  void setAK8PFHT750(const Bool_t value) {AK8PFHT750_ = value;}
  void setAK8PFHT800(const Bool_t value) {AK8PFHT800_ = value;}
  void setAK8DiPFJet300(const Bool_t value) {AK8DiPFJet300_ = value;}
  void setAK8DiPFJet280(const Bool_t value) {AK8DiPFJet280_ = value;}
  void setAK8DiPFJet300(const Bool_t value) {AK8DiPFJet300_ = value;}
  void setAK8DiPFJet300(const Bool_t value) {AK8DiPFJet300_ = value;}
  void setAK8PFHT700(const Bool_t value) {AK8PFHT700_ = value;}
  void setAK8PFHT650(const Bool_t value) {AK8PFHT650_ = value;}
  void setAK8PFHT600(const Bool_t value) {AK8PFHT600_ = value;}
  void setAK8DiPFJet280(const Bool_t value) {AK8DiPFJet280_ = value;}
  void setAK8DiPFJet250(const Bool_t value) {AK8DiPFJet250_ = value;}
  void setAK8DiPFJet280(const Bool_t value) {AK8DiPFJet280_ = value;}
  void setAK8DiPFJet250(const Bool_t value) {AK8DiPFJet250_ = value;}
  void setCaloJet260(const Bool_t value) {CaloJet260_ = value;}
  void setCaloJet500(const Bool_t value) {CaloJet500_ = value;}
  void setDimuon13(const Bool_t value) {Dimuon13_ = value;}
  void setDimuon13(const Bool_t value) {Dimuon13_ = value;}
  void setDimuon20(const Bool_t value) {Dimuon20_ = value;}
  void setDoubleEle24(const Bool_t value) {DoubleEle24_ = value;}
  void setDoubleEle25(const Bool_t value) {DoubleEle25_ = value;}
  void setDoubleEle33(const Bool_t value) {DoubleEle33_ = value;}
  void setDoubleEle33(const Bool_t value) {DoubleEle33_ = value;}
  void setDoubleEle33(const Bool_t value) {DoubleEle33_ = value;}
  void setDoubleEle33(const Bool_t value) {DoubleEle33_ = value;}
  void setDoubleMediumCombinedIsoPFTau35(const Bool_t value) {DoubleMediumCombinedIsoPFTau35_ = value;}
  void setDoubleTightCombinedIsoPFTau35(const Bool_t value) {DoubleTightCombinedIsoPFTau35_ = value;}
  void setDoubleMediumCombinedIsoPFTau40(const Bool_t value) {DoubleMediumCombinedIsoPFTau40_ = value;}
  void setDoubleTightCombinedIsoPFTau40(const Bool_t value) {DoubleTightCombinedIsoPFTau40_ = value;}
  void setDoubleMediumCombinedIsoPFTau40(const Bool_t value) {DoubleMediumCombinedIsoPFTau40_ = value;}
  void setDoubleTightCombinedIsoPFTau40(const Bool_t value) {DoubleTightCombinedIsoPFTau40_ = value;}
  void setDoubleMediumIsoPFTau35(const Bool_t value) {DoubleMediumIsoPFTau35_ = value;}
  void setDoubleMediumIsoPFTau40(const Bool_t value) {DoubleMediumIsoPFTau40_ = value;}
  void setDoubleMediumIsoPFTau40(const Bool_t value) {DoubleMediumIsoPFTau40_ = value;}
  void setDoubleEle37(const Bool_t value) {DoubleEle37_ = value;}
  void setDoubleMu33NoFiltersNoVtx(const Bool_t value) {DoubleMu33NoFiltersNoVtx_ = value;}
  void setDoubleMu38NoFiltersNoVtx(const Bool_t value) {DoubleMu38NoFiltersNoVtx_ = value;}
  void setDoubleMu23NoFiltersNoVtxDisplaced(const Bool_t value) {DoubleMu23NoFiltersNoVtxDisplaced_ = value;}
  void setDoubleMu28NoFiltersNoVtxDisplaced(const Bool_t value) {DoubleMu28NoFiltersNoVtxDisplaced_ = value;}
  void setDoubleMu0(const Bool_t value) {DoubleMu0_ = value;}
  void setDoubleMu4(const Bool_t value) {DoubleMu4_ = value;}
  void setDoubleMu4(const Bool_t value) {DoubleMu4_ = value;}
  void setDoubleMu4(const Bool_t value) {DoubleMu4_ = value;}
  void setDoubleMu4(const Bool_t value) {DoubleMu4_ = value;}
  void setDoubleMu3(const Bool_t value) {DoubleMu3_ = value;}
  void setDoubleMu4(const Bool_t value) {DoubleMu4_ = value;}
  void setMu7p5(const Bool_t value) {Mu7p5_ = value;}
  void setMu7p5(const Bool_t value) {Mu7p5_ = value;}
  void setMu7p5(const Bool_t value) {Mu7p5_ = value;}
  void setMu7p5(const Bool_t value) {Mu7p5_ = value;}
  void setMu7p5(const Bool_t value) {Mu7p5_ = value;}
  void setMu7p5(const Bool_t value) {Mu7p5_ = value;}
  void setMu7p5(const Bool_t value) {Mu7p5_ = value;}
  void setMu7p5(const Bool_t value) {Mu7p5_ = value;}
  void setDimuon0er16(const Bool_t value) {Dimuon0er16_ = value;}
  void setDimuon0er16(const Bool_t value) {Dimuon0er16_ = value;}
  void setDimuon6(const Bool_t value) {Dimuon6_ = value;}
  void setPhoton150(const Bool_t value) {Photon150_ = value;}
  void setPhoton90(const Bool_t value) {Photon90_ = value;}
  void setHT250(const Bool_t value) {HT250_ = value;}
  void setDoublePhoton60(const Bool_t value) {DoublePhoton60_ = value;}
  void setDoublePhoton85(const Bool_t value) {DoublePhoton85_ = value;}
  void setEle17(const Bool_t value) {Ele17_ = value;}
  void setEle20(const Bool_t value) {Ele20_ = value;}
  void setEle22(const Bool_t value) {Ele22_ = value;}
  void setEle22(const Bool_t value) {Ele22_ = value;}
  void setEle22(const Bool_t value) {Ele22_ = value;}
  void setEle23(const Bool_t value) {Ele23_ = value;}
  void setEle23(const Bool_t value) {Ele23_ = value;}
  void setEle24(const Bool_t value) {Ele24_ = value;}
  void setEle24(const Bool_t value) {Ele24_ = value;}
  void setEle24(const Bool_t value) {Ele24_ = value;}
  void setEle24(const Bool_t value) {Ele24_ = value;}
  void setEle25(const Bool_t value) {Ele25_ = value;}
  void setEle25(const Bool_t value) {Ele25_ = value;}
  void setEle25(const Bool_t value) {Ele25_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setEle30(const Bool_t value) {Ele30_ = value;}
  void setEle30(const Bool_t value) {Ele30_ = value;}
  void setEle30(const Bool_t value) {Ele30_ = value;}
  void setEle32(const Bool_t value) {Ele32_ = value;}
  void setEle32(const Bool_t value) {Ele32_ = value;}
  void setEle32(const Bool_t value) {Ele32_ = value;}
  void setEle32(const Bool_t value) {Ele32_ = value;}
  void setEle35(const Bool_t value) {Ele35_ = value;}
  void setEle35(const Bool_t value) {Ele35_ = value;}
  void setEle36(const Bool_t value) {Ele36_ = value;}
  void setEle45(const Bool_t value) {Ele45_ = value;}
  void setEle45(const Bool_t value) {Ele45_ = value;}
  void setEle45(const Bool_t value) {Ele45_ = value;}
  void setEle105(const Bool_t value) {Ele105_ = value;}
  void setEle30WP60(const Bool_t value) {Ele30WP60_ = value;}
  void setEle30WP60(const Bool_t value) {Ele30WP60_ = value;}
  void setHT200(const Bool_t value) {HT200_ = value;}
  void setHT275(const Bool_t value) {HT275_ = value;}
  void setHT325(const Bool_t value) {HT325_ = value;}
  void setHT425(const Bool_t value) {HT425_ = value;}
  void setHT575(const Bool_t value) {HT575_ = value;}
  void setHT410to430(const Bool_t value) {HT410to430_ = value;}
  void setHT430to450(const Bool_t value) {HT430to450_ = value;}
  void setHT450to470(const Bool_t value) {HT450to470_ = value;}
  void setHT470to500(const Bool_t value) {HT470to500_ = value;}
  void setHT500to550(const Bool_t value) {HT500to550_ = value;}
  void setHT550to650(const Bool_t value) {HT550to650_ = value;}
  void setHT650(const Bool_t value) {HT650_ = value;}
  void setMu16(const Bool_t value) {Mu16_ = value;}
  void setIsoMu16(const Bool_t value) {IsoMu16_ = value;}
  void setIsoMu16(const Bool_t value) {IsoMu16_ = value;}
  void setIsoMu17(const Bool_t value) {IsoMu17_ = value;}
  void setIsoMu17(const Bool_t value) {IsoMu17_ = value;}
  void setIsoMu17(const Bool_t value) {IsoMu17_ = value;}
  void setDoubleIsoMu17(const Bool_t value) {DoubleIsoMu17_ = value;}
  void setDoubleIsoMu17(const Bool_t value) {DoubleIsoMu17_ = value;}
  void setIsoMu18(const Bool_t value) {IsoMu18_ = value;}
  void setIsoMu19(const Bool_t value) {IsoMu19_ = value;}
  void setIsoMu19(const Bool_t value) {IsoMu19_ = value;}
  void setIsoMu19(const Bool_t value) {IsoMu19_ = value;}
  void setIsoMu19(const Bool_t value) {IsoMu19_ = value;}
  void setIsoMu19(const Bool_t value) {IsoMu19_ = value;}
  void setIsoMu19(const Bool_t value) {IsoMu19_ = value;}
  void setIsoMu21(const Bool_t value) {IsoMu21_ = value;}
  void setIsoMu21(const Bool_t value) {IsoMu21_ = value;}
  void setIsoMu20(const Bool_t value) {IsoMu20_ = value;}
  void setIsoMu21(const Bool_t value) {IsoMu21_ = value;}
  void setIsoMu21(const Bool_t value) {IsoMu21_ = value;}
  void setIsoMu21(const Bool_t value) {IsoMu21_ = value;}
  void setIsoMu22(const Bool_t value) {IsoMu22_ = value;}
  void setIsoMu22(const Bool_t value) {IsoMu22_ = value;}
  void setIsoMu24(const Bool_t value) {IsoMu24_ = value;}
  void setIsoMu27(const Bool_t value) {IsoMu27_ = value;}
  void setIsoTkMu18(const Bool_t value) {IsoTkMu18_ = value;}
  void setIsoTkMu20(const Bool_t value) {IsoTkMu20_ = value;}
  void setIsoTkMu22(const Bool_t value) {IsoTkMu22_ = value;}
  void setIsoTkMu22(const Bool_t value) {IsoTkMu22_ = value;}
  void setIsoTkMu24(const Bool_t value) {IsoTkMu24_ = value;}
  void setIsoTkMu27(const Bool_t value) {IsoTkMu27_ = value;}
  void setJetE30(const Bool_t value) {JetE30_ = value;}
  void setJetE30(const Bool_t value) {JetE30_ = value;}
  void setJetE50(const Bool_t value) {JetE50_ = value;}
  void setJetE70(const Bool_t value) {JetE70_ = value;}
  void setL1SingleMu18(const Bool_t value) {L1SingleMu18_ = value;}
  void setL2Mu10(const Bool_t value) {L2Mu10_ = value;}
  void setL1SingleMuOpen(const Bool_t value) {L1SingleMuOpen_ = value;}
  void setL1SingleMuOpen(const Bool_t value) {L1SingleMuOpen_ = value;}
  void setL2DoubleMu23(const Bool_t value) {L2DoubleMu23_ = value;}
  void setL2DoubleMu28(const Bool_t value) {L2DoubleMu28_ = value;}
  void setL2DoubleMu38(const Bool_t value) {L2DoubleMu38_ = value;}
  void setL2Mu10(const Bool_t value) {L2Mu10_ = value;}
  void setL2Mu10(const Bool_t value) {L2Mu10_ = value;}
  void setL2Mu45(const Bool_t value) {L2Mu45_ = value;}
  void setL2Mu40(const Bool_t value) {L2Mu40_ = value;}
  void setLooseIsoPFTau50(const Bool_t value) {LooseIsoPFTau50_ = value;}
  void setLooseIsoPFTau50(const Bool_t value) {LooseIsoPFTau50_ = value;}
  void setLooseIsoPFTau50(const Bool_t value) {LooseIsoPFTau50_ = value;}
  void setLooseIsoPFTau50(const Bool_t value) {LooseIsoPFTau50_ = value;}
  void setLooseIsoPFTau50(const Bool_t value) {LooseIsoPFTau50_ = value;}
  void setPFTau120(const Bool_t value) {PFTau120_ = value;}
  void setPFTau140(const Bool_t value) {PFTau140_ = value;}
  void setVLooseIsoPFTau120(const Bool_t value) {VLooseIsoPFTau120_ = value;}
  void setVLooseIsoPFTau140(const Bool_t value) {VLooseIsoPFTau140_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu20(const Bool_t value) {Mu20_ = value;}
  void setMu20(const Bool_t value) {Mu20_ = value;}
  void setMu20(const Bool_t value) {Mu20_ = value;}
  void setMu20(const Bool_t value) {Mu20_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu25(const Bool_t value) {Mu25_ = value;}
  void setMu27(const Bool_t value) {Mu27_ = value;}
  void setMu30(const Bool_t value) {Mu30_ = value;}
  void setMu30(const Bool_t value) {Mu30_ = value;}
  void setMu40(const Bool_t value) {Mu40_ = value;}
  void setMu40(const Bool_t value) {Mu40_ = value;}
  void setMu20(const Bool_t value) {Mu20_ = value;}
  void setTkMu17(const Bool_t value) {TkMu17_ = value;}
  void setTkMu17(const Bool_t value) {TkMu17_ = value;}
  void setTkMu17(const Bool_t value) {TkMu17_ = value;}
  void setTkMu20(const Bool_t value) {TkMu20_ = value;}
  void setMu24(const Bool_t value) {Mu24_ = value;}
  void setTkMu24(const Bool_t value) {TkMu24_ = value;}
  void setMu27(const Bool_t value) {Mu27_ = value;}
  void setTkMu27(const Bool_t value) {TkMu27_ = value;}
  void setMu45(const Bool_t value) {Mu45_ = value;}
  void setMu50(const Bool_t value) {Mu50_ = value;}
  void setTkMu50(const Bool_t value) {TkMu50_ = value;}
  void setMu38NoFiltersNoVtx(const Bool_t value) {Mu38NoFiltersNoVtx_ = value;}
  void setMu42NoFiltersNoVtx(const Bool_t value) {Mu42NoFiltersNoVtx_ = value;}
  void setMu28NoFiltersNoVtxDisplaced(const Bool_t value) {Mu28NoFiltersNoVtxDisplaced_ = value;}
  void setMu33NoFiltersNoVtxDisplaced(const Bool_t value) {Mu33NoFiltersNoVtxDisplaced_ = value;}
  void setMu23NoFiltersNoVtx(const Bool_t value) {Mu23NoFiltersNoVtx_ = value;}
  void setDoubleMu18NoFiltersNoVtx(const Bool_t value) {DoubleMu18NoFiltersNoVtx_ = value;}
  void setMu33NoFiltersNoVtxDisplaced(const Bool_t value) {Mu33NoFiltersNoVtxDisplaced_ = value;}
  void setMu33NoFiltersNoVtxDisplaced(const Bool_t value) {Mu33NoFiltersNoVtxDisplaced_ = value;}
  void setMu28NoFiltersNoVtx(const Bool_t value) {Mu28NoFiltersNoVtx_ = value;}
  void setMu38NoFiltersNoVtxDisplaced(const Bool_t value) {Mu38NoFiltersNoVtxDisplaced_ = value;}
  void setMu38NoFiltersNoVtxDisplaced(const Bool_t value) {Mu38NoFiltersNoVtxDisplaced_ = value;}
  void setMu38NoFiltersNoVtx(const Bool_t value) {Mu38NoFiltersNoVtx_ = value;}
  void setMu28NoFiltersNoVtx(const Bool_t value) {Mu28NoFiltersNoVtx_ = value;}
  void setPFHT300(const Bool_t value) {PFHT300_ = value;}
  void setPFHT300(const Bool_t value) {PFHT300_ = value;}
  void setPFHT550(const Bool_t value) {PFHT550_ = value;}
  void setPFHT650(const Bool_t value) {PFHT650_ = value;}
  void setPFHT750(const Bool_t value) {PFHT750_ = value;}
  void setPFHT750(const Bool_t value) {PFHT750_ = value;}
  void setPFHT750(const Bool_t value) {PFHT750_ = value;}
  void setPFHT800(const Bool_t value) {PFHT800_ = value;}
  void setPFHT850(const Bool_t value) {PFHT850_ = value;}
  void setPFJet15(const Bool_t value) {PFJet15_ = value;}
  void setPFJet25(const Bool_t value) {PFJet25_ = value;}
  void setDiPFJet15(const Bool_t value) {DiPFJet15_ = value;}
  void setDiPFJet25(const Bool_t value) {DiPFJet25_ = value;}
  void setDiPFJet15(const Bool_t value) {DiPFJet15_ = value;}
  void setDiPFJet25(const Bool_t value) {DiPFJet25_ = value;}
  void setDiPFJetAve15(const Bool_t value) {DiPFJetAve15_ = value;}
  void setDiPFJetAve25(const Bool_t value) {DiPFJetAve25_ = value;}
  void setDiPFJetAve35(const Bool_t value) {DiPFJetAve35_ = value;}
  void setAK8PFJet40(const Bool_t value) {AK8PFJet40_ = value;}
  void setAK8PFJet60(const Bool_t value) {AK8PFJet60_ = value;}
  void setAK8PFJet80(const Bool_t value) {AK8PFJet80_ = value;}
  void setAK8PFJet140(const Bool_t value) {AK8PFJet140_ = value;}
  void setAK8PFJet200(const Bool_t value) {AK8PFJet200_ = value;}
  void setAK8PFJet260(const Bool_t value) {AK8PFJet260_ = value;}
  void setAK8PFJet320(const Bool_t value) {AK8PFJet320_ = value;}
  void setAK8PFJet400(const Bool_t value) {AK8PFJet400_ = value;}
  void setAK8PFJet450(const Bool_t value) {AK8PFJet450_ = value;}
  void setAK8PFJet500(const Bool_t value) {AK8PFJet500_ = value;}
  void setPFJet40(const Bool_t value) {PFJet40_ = value;}
  void setPFJet60(const Bool_t value) {PFJet60_ = value;}
  void setPFJet80(const Bool_t value) {PFJet80_ = value;}
  void setPFJet140(const Bool_t value) {PFJet140_ = value;}
  void setPFJet200(const Bool_t value) {PFJet200_ = value;}
  void setPFJet260(const Bool_t value) {PFJet260_ = value;}
  void setPFJet320(const Bool_t value) {PFJet320_ = value;}
  void setPFJet400(const Bool_t value) {PFJet400_ = value;}
  void setPFJet450(const Bool_t value) {PFJet450_ = value;}
  void setPFJet500(const Bool_t value) {PFJet500_ = value;}
  void setDiPFJetAve40(const Bool_t value) {DiPFJetAve40_ = value;}
  void setDiPFJetAve60(const Bool_t value) {DiPFJetAve60_ = value;}
  void setDiPFJetAve80(const Bool_t value) {DiPFJetAve80_ = value;}
  void setDiPFJetAve140(const Bool_t value) {DiPFJetAve140_ = value;}
  void setDiPFJetAve200(const Bool_t value) {DiPFJetAve200_ = value;}
  void setDiPFJetAve260(const Bool_t value) {DiPFJetAve260_ = value;}
  void setDiPFJetAve320(const Bool_t value) {DiPFJetAve320_ = value;}
  void setDiPFJetAve400(const Bool_t value) {DiPFJetAve400_ = value;}
  void setDiPFJetAve500(const Bool_t value) {DiPFJetAve500_ = value;}
  void setDiPFJetAve60(const Bool_t value) {DiPFJetAve60_ = value;}
  void setDiPFJetAve80(const Bool_t value) {DiPFJetAve80_ = value;}
  void setDiPFJetAve100(const Bool_t value) {DiPFJetAve100_ = value;}
  void setDiPFJetAve160(const Bool_t value) {DiPFJetAve160_ = value;}
  void setDiPFJetAve220(const Bool_t value) {DiPFJetAve220_ = value;}
  void setDiPFJetAve300(const Bool_t value) {DiPFJetAve300_ = value;}
  void setDiPFJet40(const Bool_t value) {DiPFJet40_ = value;}
  void setDiPFJet40(const Bool_t value) {DiPFJet40_ = value;}
  void setDiCentralPFJet170(const Bool_t value) {DiCentralPFJet170_ = value;}
  void setSingleCentralPFJet170(const Bool_t value) {SingleCentralPFJet170_ = value;}
  void setDiCentralPFJet170(const Bool_t value) {DiCentralPFJet170_ = value;}
  void setDiCentralPFJet220(const Bool_t value) {DiCentralPFJet220_ = value;}
  void setDiCentralPFJet330(const Bool_t value) {DiCentralPFJet330_ = value;}
  void setDiCentralPFJet430(const Bool_t value) {DiCentralPFJet430_ = value;}
  void setPFHT125(const Bool_t value) {PFHT125_ = value;}
  void setPFHT200(const Bool_t value) {PFHT200_ = value;}
  void setPFHT250(const Bool_t value) {PFHT250_ = value;}
  void setPFHT300(const Bool_t value) {PFHT300_ = value;}
  void setPFHT350(const Bool_t value) {PFHT350_ = value;}
  void setPFHT400(const Bool_t value) {PFHT400_ = value;}
  void setPFHT475(const Bool_t value) {PFHT475_ = value;}
  void setPFHT600(const Bool_t value) {PFHT600_ = value;}
  void setPFHT650(const Bool_t value) {PFHT650_ = value;}
  void setPFHT800(const Bool_t value) {PFHT800_ = value;}
  void setPFHT900(const Bool_t value) {PFHT900_ = value;}
  void setPFHT200(const Bool_t value) {PFHT200_ = value;}
  void setPFHT200(const Bool_t value) {PFHT200_ = value;}
  void setPFHT200(const Bool_t value) {PFHT200_ = value;}
  void setPFHT250(const Bool_t value) {PFHT250_ = value;}
  void setPFHT250(const Bool_t value) {PFHT250_ = value;}
  void setPFHT300(const Bool_t value) {PFHT300_ = value;}
  void setPFHT300(const Bool_t value) {PFHT300_ = value;}
  void setPFHT350(const Bool_t value) {PFHT350_ = value;}
  void setPFHT350(const Bool_t value) {PFHT350_ = value;}
  void setPFHT400(const Bool_t value) {PFHT400_ = value;}
  void setPFHT400(const Bool_t value) {PFHT400_ = value;}
  void setMET60(const Bool_t value) {MET60_ = value;}
  void setMET75(const Bool_t value) {MET75_ = value;}
  void setMET90(const Bool_t value) {MET90_ = value;}
  void setPFMET120(const Bool_t value) {PFMET120_ = value;}
  void setPFMET120(const Bool_t value) {PFMET120_ = value;}
  void setPFMET170(const Bool_t value) {PFMET170_ = value;}
  void setPFMET170(const Bool_t value) {PFMET170_ = value;}
  void setPFMET170(const Bool_t value) {PFMET170_ = value;}
  void setPFMET170(const Bool_t value) {PFMET170_ = value;}
  void setPFMET170(const Bool_t value) {PFMET170_ = value;}
  void setPFMET170(const Bool_t value) {PFMET170_ = value;}
  void setPFMETTypeOne190(const Bool_t value) {PFMETTypeOne190_ = value;}
  void setPFMET90(const Bool_t value) {PFMET90_ = value;}
  void setPFMET100(const Bool_t value) {PFMET100_ = value;}
  void setPFMET100(const Bool_t value) {PFMET100_ = value;}
  void setPFMET110(const Bool_t value) {PFMET110_ = value;}
  void setPFMET120(const Bool_t value) {PFMET120_ = value;}
  void setCaloMHTNoPU90(const Bool_t value) {CaloMHTNoPU90_ = value;}
  void setCaloMHTNoPU90(const Bool_t value) {CaloMHTNoPU90_ = value;}
  void setQuadPFJet(const Bool_t value) {QuadPFJet_ = value;}
  void setQuadPFJet(const Bool_t value) {QuadPFJet_ = value;}
  void setQuadPFJet(const Bool_t value) {QuadPFJet_ = value;}
  void setQuadPFJet(const Bool_t value) {QuadPFJet_ = value;}
  void setQuadPFJet(const Bool_t value) {QuadPFJet_ = value;}
  void setL1(const Bool_t value) {L1_ = value;}
  void setQuadJet45(const Bool_t value) {QuadJet45_ = value;}
  void setQuadJet45(const Bool_t value) {QuadJet45_ = value;}
  void setDoubleJet90(const Bool_t value) {DoubleJet90_ = value;}
  void setDoubleJet90(const Bool_t value) {DoubleJet90_ = value;}
  void setDoubleJetsC100(const Bool_t value) {DoubleJetsC100_ = value;}
  void setDoubleJetsC100(const Bool_t value) {DoubleJetsC100_ = value;}
  void setDoubleJetsC112(const Bool_t value) {DoubleJetsC112_ = value;}
  void setDoubleJetsC112(const Bool_t value) {DoubleJetsC112_ = value;}
  void setDoubleJetsC100(const Bool_t value) {DoubleJetsC100_ = value;}
  void setDoubleJetsC100(const Bool_t value) {DoubleJetsC100_ = value;}
  void setDoubleJetsC100(const Bool_t value) {DoubleJetsC100_ = value;}
  void setDoubleJetsC100(const Bool_t value) {DoubleJetsC100_ = value;}
  void setPhoton135(const Bool_t value) {Photon135_ = value;}
  void setPhoton20(const Bool_t value) {Photon20_ = value;}
  void setPhoton22(const Bool_t value) {Photon22_ = value;}
  void setPhoton22(const Bool_t value) {Photon22_ = value;}
  void setPhoton250(const Bool_t value) {Photon250_ = value;}
  void setPhoton300(const Bool_t value) {Photon300_ = value;}
  void setPhoton26(const Bool_t value) {Photon26_ = value;}
  void setPhoton36(const Bool_t value) {Photon36_ = value;}
  void setPhoton36(const Bool_t value) {Photon36_ = value;}
  void setPhoton36(const Bool_t value) {Photon36_ = value;}
  void setPhoton50(const Bool_t value) {Photon50_ = value;}
  void setPhoton50(const Bool_t value) {Photon50_ = value;}
  void setPhoton75(const Bool_t value) {Photon75_ = value;}
  void setPhoton75(const Bool_t value) {Photon75_ = value;}
  void setPhoton90(const Bool_t value) {Photon90_ = value;}
  void setPhoton90(const Bool_t value) {Photon90_ = value;}
  void setPhoton120(const Bool_t value) {Photon120_ = value;}
  void setPhoton120(const Bool_t value) {Photon120_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setEle8(const Bool_t value) {Ele8_ = value;}
  void setEle12(const Bool_t value) {Ele12_ = value;}
  void setEle17(const Bool_t value) {Ele17_ = value;}
  void setEle23(const Bool_t value) {Ele23_ = value;}
  void setBTagMu(const Bool_t value) {BTagMu_ = value;}
  void setBTagMu(const Bool_t value) {BTagMu_ = value;}
  void setBTagMu(const Bool_t value) {BTagMu_ = value;}
  void setBTagMu(const Bool_t value) {BTagMu_ = value;}
  void setBTagMu(const Bool_t value) {BTagMu_ = value;}
  void setBTagMu(const Bool_t value) {BTagMu_ = value;}
  void setBTagMu(const Bool_t value) {BTagMu_ = value;}
  void setEle23(const Bool_t value) {Ele23_ = value;}
  void setEle23(const Bool_t value) {Ele23_ = value;}
  void setEle17(const Bool_t value) {Ele17_ = value;}
  void setEle16(const Bool_t value) {Ele16_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setMu12(const Bool_t value) {Mu12_ = value;}
  void setMu12(const Bool_t value) {Mu12_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu23(const Bool_t value) {Mu23_ = value;}
  void setMu23(const Bool_t value) {Mu23_ = value;}
  void setMu23(const Bool_t value) {Mu23_ = value;}
  void setMu23(const Bool_t value) {Mu23_ = value;}
  void setMu30(const Bool_t value) {Mu30_ = value;}
  void setMu33(const Bool_t value) {Mu33_ = value;}
  void setMu37(const Bool_t value) {Mu37_ = value;}
  void setMu27(const Bool_t value) {Mu27_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setMu12(const Bool_t value) {Mu12_ = value;}
  void setMu12(const Bool_t value) {Mu12_ = value;}
  void setMu12(const Bool_t value) {Mu12_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setDiMu9(const Bool_t value) {DiMu9_ = value;}
  void setTripleMu(const Bool_t value) {TripleMu_ = value;}
  void setTripleMu(const Bool_t value) {TripleMu_ = value;}
  void setMu3er(const Bool_t value) {Mu3er_ = value;}
  void setMu6(const Bool_t value) {Mu6_ = value;}
  void setMu6(const Bool_t value) {Mu6_ = value;}
  void setMu14er(const Bool_t value) {Mu14er_ = value;}
  void setEle17(const Bool_t value) {Ele17_ = value;}
  void setEle23(const Bool_t value) {Ele23_ = value;}
  void setEle12(const Bool_t value) {Ele12_ = value;}
  void setEle17(const Bool_t value) {Ele17_ = value;}
  void setEle17(const Bool_t value) {Ele17_ = value;}
  void setEle23(const Bool_t value) {Ele23_ = value;}
  void setPFHT650(const Bool_t value) {PFHT650_ = value;}
  void setPFHT650(const Bool_t value) {PFHT650_ = value;}
  void setPhoton22(const Bool_t value) {Photon22_ = value;}
  void setPhoton30(const Bool_t value) {Photon30_ = value;}
  void setPhoton36(const Bool_t value) {Photon36_ = value;}
  void setPhoton50(const Bool_t value) {Photon50_ = value;}
  void setPhoton75(const Bool_t value) {Photon75_ = value;}
  void setPhoton90(const Bool_t value) {Photon90_ = value;}
  void setPhoton120(const Bool_t value) {Photon120_ = value;}
  void setPhoton175(const Bool_t value) {Photon175_ = value;}
  void setPhoton165(const Bool_t value) {Photon165_ = value;}
  void setPhoton22(const Bool_t value) {Photon22_ = value;}
  void setPhoton30(const Bool_t value) {Photon30_ = value;}
  void setPhoton36(const Bool_t value) {Photon36_ = value;}
  void setPhoton50(const Bool_t value) {Photon50_ = value;}
  void setPhoton75(const Bool_t value) {Photon75_ = value;}
  void setPhoton90(const Bool_t value) {Photon90_ = value;}
  void setPhoton120(const Bool_t value) {Photon120_ = value;}
  void setPhoton165(const Bool_t value) {Photon165_ = value;}
  void setDiphoton30(const Bool_t value) {Diphoton30_ = value;}
  void setDiphoton30(const Bool_t value) {Diphoton30_ = value;}
  void setDiphoton30PV(const Bool_t value) {Diphoton30PV_ = value;}
  void setDiphoton30(const Bool_t value) {Diphoton30_ = value;}
  void setDiphoton30EB(const Bool_t value) {Diphoton30EB_ = value;}
  void setDimuon0(const Bool_t value) {Dimuon0_ = value;}
  void setDimuon0(const Bool_t value) {Dimuon0_ = value;}
  void setQuadMuon0(const Bool_t value) {QuadMuon0_ = value;}
  void setQuadMuon0(const Bool_t value) {QuadMuon0_ = value;}
  void setRsq0p25(const Bool_t value) {Rsq0p25_ = value;}
  void setRsqMR240(const Bool_t value) {RsqMR240_ = value;}
  void setRsqMR240(const Bool_t value) {RsqMR240_ = value;}
  void setRsq0p25(const Bool_t value) {Rsq0p25_ = value;}
  void setRsq0p30(const Bool_t value) {Rsq0p30_ = value;}
  void setRsqMR240(const Bool_t value) {RsqMR240_ = value;}
  void setRsqMR240(const Bool_t value) {RsqMR240_ = value;}
  void setRsqMR270(const Bool_t value) {RsqMR270_ = value;}
  void setRsqMR270(const Bool_t value) {RsqMR270_ = value;}
  void setRsq0p02(const Bool_t value) {Rsq0p02_ = value;}
  void setRsq0p02(const Bool_t value) {Rsq0p02_ = value;}
  void setRsq0p02(const Bool_t value) {Rsq0p02_ = value;}
  void setRsq0p02(const Bool_t value) {Rsq0p02_ = value;}
  void setRsq0p02(const Bool_t value) {Rsq0p02_ = value;}
  void setHT200(const Bool_t value) {HT200_ = value;}
  void setHT250(const Bool_t value) {HT250_ = value;}
  void setHT350(const Bool_t value) {HT350_ = value;}
  void setHT350(const Bool_t value) {HT350_ = value;}
  void setHT350(const Bool_t value) {HT350_ = value;}
  void setHT350(const Bool_t value) {HT350_ = value;}
  void setHT400(const Bool_t value) {HT400_ = value;}
  void setHT500(const Bool_t value) {HT500_ = value;}
  void setHT550(const Bool_t value) {HT550_ = value;}
  void setHT550(const Bool_t value) {HT550_ = value;}
  void setHT650(const Bool_t value) {HT650_ = value;}
  void setHT750(const Bool_t value) {HT750_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setVBF(const Bool_t value) {VBF_ = value;}
  void setPFMETNoMu90(const Bool_t value) {PFMETNoMu90_ = value;}
  void setPFMETNoMu100(const Bool_t value) {PFMETNoMu100_ = value;}
  void setPFMETNoMu110(const Bool_t value) {PFMETNoMu110_ = value;}
  void setPFMETNoMu120(const Bool_t value) {PFMETNoMu120_ = value;}
  void setMonoCentralPFJet80(const Bool_t value) {MonoCentralPFJet80_ = value;}
  void setMonoCentralPFJet80(const Bool_t value) {MonoCentralPFJet80_ = value;}
  void setMonoCentralPFJet80(const Bool_t value) {MonoCentralPFJet80_ = value;}
  void setMonoCentralPFJet80(const Bool_t value) {MonoCentralPFJet80_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setPhoton90(const Bool_t value) {Photon90_ = value;}
  void setDoubleMu8(const Bool_t value) {DoubleMu8_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setDoubleEle8(const Bool_t value) {DoubleEle8_ = value;}
  void setDoubleMu8(const Bool_t value) {DoubleMu8_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setDoubleEle8(const Bool_t value) {DoubleEle8_ = value;}
  void setMu10(const Bool_t value) {Mu10_ = value;}
  void setDoubleMu3(const Bool_t value) {DoubleMu3_ = value;}
  void setEle10(const Bool_t value) {Ele10_ = value;}
  void setEle15(const Bool_t value) {Ele15_ = value;}
  void setEle15(const Bool_t value) {Ele15_ = value;}
  void setEle15(const Bool_t value) {Ele15_ = value;}
  void setEle15(const Bool_t value) {Ele15_ = value;}
  void setEle15(const Bool_t value) {Ele15_ = value;}
  void setEle15(const Bool_t value) {Ele15_ = value;}
  void setEle50(const Bool_t value) {Ele50_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setMu10(const Bool_t value) {Mu10_ = value;}
  void setMu15(const Bool_t value) {Mu15_ = value;}
  void setMu15(const Bool_t value) {Mu15_ = value;}
  void setMu15(const Bool_t value) {Mu15_ = value;}
  void setMu15(const Bool_t value) {Mu15_ = value;}
  void setMu15(const Bool_t value) {Mu15_ = value;}
  void setMu15(const Bool_t value) {Mu15_ = value;}
  void setMu50(const Bool_t value) {Mu50_ = value;}
  void setDimuon16(const Bool_t value) {Dimuon16_ = value;}
  void setDimuon10(const Bool_t value) {Dimuon10_ = value;}
  void setDimuon8(const Bool_t value) {Dimuon8_ = value;}
  void setDimuon8(const Bool_t value) {Dimuon8_ = value;}
  void setDimuon0(const Bool_t value) {Dimuon0_ = value;}
  void setMu16(const Bool_t value) {Mu16_ = value;}
  void setMu16(const Bool_t value) {Mu16_ = value;}
  void setTrkMu15(const Bool_t value) {TrkMu15_ = value;}
  void setTrkMu17(const Bool_t value) {TrkMu17_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu3(const Bool_t value) {Mu3_ = value;}
  void setEle8(const Bool_t value) {Ele8_ = value;}
  void setEle12(const Bool_t value) {Ele12_ = value;}
  void setEle17(const Bool_t value) {Ele17_ = value;}
  void setEle23(const Bool_t value) {Ele23_ = value;}
  void setEle50(const Bool_t value) {Ele50_ = value;}
  void setEle50(const Bool_t value) {Ele50_ = value;}
  void setPFHT400(const Bool_t value) {PFHT400_ = value;}
  void setPFHT450(const Bool_t value) {PFHT450_ = value;}
  void setPFHT400(const Bool_t value) {PFHT400_ = value;}
  void setPFHT450(const Bool_t value) {PFHT450_ = value;}
  void setEle115(const Bool_t value) {Ele115_ = value;}
  void setMu55(const Bool_t value) {Mu55_ = value;}
  void setPhoton42(const Bool_t value) {Photon42_ = value;}
  void setPhoton90(const Bool_t value) {Photon90_ = value;}
  void setPixelTracks(const Bool_t value) {PixelTracks_ = value;}
  void setPixelTracks(const Bool_t value) {PixelTracks_ = value;}
  void setPixelTracks(const Bool_t value) {PixelTracks_ = value;}
  void setPixelTracks(const Bool_t value) {PixelTracks_ = value;}
  void setPixelTracks(const Bool_t value) {PixelTracks_ = value;}
  void setFullTracks(const Bool_t value) {FullTracks_ = value;}
  void setFullTracks(const Bool_t value) {FullTracks_ = value;}
  void setFullTracks(const Bool_t value) {FullTracks_ = value;}
  void setFullTracks(const Bool_t value) {FullTracks_ = value;}
  void setECALHT800(const Bool_t value) {ECALHT800_ = value;}
  void setDiSC30(const Bool_t value) {DiSC30_ = value;}
  void setPhoton125(const Bool_t value) {Photon125_ = value;}
  void setMET100(const Bool_t value) {MET100_ = value;}
  void setMET150(const Bool_t value) {MET150_ = value;}
  void setMET200(const Bool_t value) {MET200_ = value;}
  void setEle27(const Bool_t value) {Ele27_ = value;}
  void setL1FatEvents(const Bool_t value) {L1FatEvents_ = value;}
  void setPhysics(const Bool_t value) {Physics_ = value;}
  void setL1FatEvents(const Bool_t value) {L1FatEvents_ = value;}
  void setL1FatEvents(const Bool_t value) {L1FatEvents_ = value;}
  void setL1FatEvents(const Bool_t value) {L1FatEvents_ = value;}
  void setL1FatEvents(const Bool_t value) {L1FatEvents_ = value;}
  void setRandom(const Bool_t value) {Random_ = value;}
  void setZeroBias(const Bool_t value) {ZeroBias_ = value;}
  void setAK4CaloJet30(const Bool_t value) {AK4CaloJet30_ = value;}
  void setAK4CaloJet40(const Bool_t value) {AK4CaloJet40_ = value;}
  void setAK4CaloJet50(const Bool_t value) {AK4CaloJet50_ = value;}
  void setAK4CaloJet80(const Bool_t value) {AK4CaloJet80_ = value;}
  void setAK4CaloJet100(const Bool_t value) {AK4CaloJet100_ = value;}
  void setAK4PFJet30(const Bool_t value) {AK4PFJet30_ = value;}
  void setAK4PFJet50(const Bool_t value) {AK4PFJet50_ = value;}
  void setAK4PFJet80(const Bool_t value) {AK4PFJet80_ = value;}
  void setAK4PFJet100(const Bool_t value) {AK4PFJet100_ = value;}
  void setHISinglePhoton10(const Bool_t value) {HISinglePhoton10_ = value;}
  void setHISinglePhoton15(const Bool_t value) {HISinglePhoton15_ = value;}
  void setHISinglePhoton20(const Bool_t value) {HISinglePhoton20_ = value;}
  void setHISinglePhoton40(const Bool_t value) {HISinglePhoton40_ = value;}
  void setHISinglePhoton60(const Bool_t value) {HISinglePhoton60_ = value;}
  void setEcalCalibration(const Bool_t value) {EcalCalibration_ = value;}
  void setHcalCalibration(const Bool_t value) {HcalCalibration_ = value;}
  void setGlobalRunHPDNoise(const Bool_t value) {GlobalRunHPDNoise_ = value;}
  void setL1BptxMinus(const Bool_t value) {L1BptxMinus_ = value;}
  void setL1BptxPlus(const Bool_t value) {L1BptxPlus_ = value;}
  void setL1NotBptxOR(const Bool_t value) {L1NotBptxOR_ = value;}
  void setL1BeamGasMinus(const Bool_t value) {L1BeamGasMinus_ = value;}
  void setL1BeamGasPlus(const Bool_t value) {L1BeamGasPlus_ = value;}
  void setL1BptxXOR(const Bool_t value) {L1BptxXOR_ = value;}
  void setL1MinimumBiasHF(const Bool_t value) {L1MinimumBiasHF_ = value;}
  void setL1MinimumBiasHF(const Bool_t value) {L1MinimumBiasHF_ = value;}
  void setHcalNZS(const Bool_t value) {HcalNZS_ = value;}
  void setHcalPhiSym(const Bool_t value) {HcalPhiSym_ = value;}
  void setHcalIsolatedbunch(const Bool_t value) {HcalIsolatedbunch_ = value;}
  void setZeroBias(const Bool_t value) {ZeroBias_ = value;}
  void setZeroBias(const Bool_t value) {ZeroBias_ = value;}
  void setZeroBias(const Bool_t value) {ZeroBias_ = value;}
  void setZeroBias(const Bool_t value) {ZeroBias_ = value;}
  void setZeroBias(const Bool_t value) {ZeroBias_ = value;}
  void setZeroBias(const Bool_t value) {ZeroBias_ = value;}
  void setPhoton500(const Bool_t value) {Photon500_ = value;}
  void setPhoton600(const Bool_t value) {Photon600_ = value;}
  void setMu300(const Bool_t value) {Mu300_ = value;}
  void setMu350(const Bool_t value) {Mu350_ = value;}
  void setMET250(const Bool_t value) {MET250_ = value;}
  void setMET300(const Bool_t value) {MET300_ = value;}
  void setMET600(const Bool_t value) {MET600_ = value;}
  void setMET700(const Bool_t value) {MET700_ = value;}
  void setPFMET300(const Bool_t value) {PFMET300_ = value;}
  void setPFMET400(const Bool_t value) {PFMET400_ = value;}
  void setPFMET500(const Bool_t value) {PFMET500_ = value;}
  void setPFMET600(const Bool_t value) {PFMET600_ = value;}
  void setEle250(const Bool_t value) {Ele250_ = value;}
  void setEle300(const Bool_t value) {Ele300_ = value;}
  void setHT2000(const Bool_t value) {HT2000_ = value;}
  void setHT2500(const Bool_t value) {HT2500_ = value;}
  void setIsoTrackHE(const Bool_t value) {IsoTrackHE_ = value;}
  void setIsoTrackHB(const Bool_t value) {IsoTrackHB_ = value;}
  void setHLTriggerFinalPath(const Bool_t value) {HLTriggerFinalPath_ = value;}
};

class S: public TLorentzVector{
friend class URStreamer;
public:
//  S(const Float_ &i_dlen_,const Float_ &i_dlenSig_,const Float_ &i_pAngle_,const Float_ &i_chi2_,const Float_ &i_ndof_,const Float_ &i_x_,const Float_ &i_y_,const Float_ &i_z_):
//    
//  {}
  S():
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
  Float_ dlen() const {return dlen_;}
  Float_ dlenSig() const {return dlenSig_;}
  Float_ pAngle() const {return pAngle_;}
  Float_ chi2() const {return chi2_;}
  Float_ ndof() const {return ndof_;}
  Float_ x() const {return x_;}
  Float_ y() const {return y_;}
  Float_ z() const {return z_;}
private:
  Float_ dlen_;
  Float_ dlenSig_;
  Float_ pAngle_;
  Float_ chi2_;
  Float_ ndof_;
  Float_ x_;
  Float_ y_;
  Float_ z_;
  void setdlen(const Float_ value) {dlen_ = value;}
  void setdlenSig(const Float_ value) {dlenSig_ = value;}
  void setpAngle(const Float_ value) {pAngle_ = value;}
  void setchi2(const Float_ value) {chi2_ = value;}
  void setndof(const Float_ value) {ndof_ = value;}
  void setx(const Float_ value) {x_ = value;}
  void sety(const Float_ value) {y_ = value;}
  void setz(const Float_ value) {z_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Met: public TLorentzVector{
friend class URStreamer;
public:
//  Met(const Float_t &i_MetUnclustEnUpDeltaX_,const Float_t &i_MetUnclustEnUpDeltaY_,const Float_t &i_covXX_,const Float_t &i_covXY_,const Float_t &i_covYY_,const Float_t &i_significance_,const Float_t &i_sumEt_,const Float_t &i_fiducialGenPhi_,const Float_t &i_fiducialGenPt_):
//    
//  {}
  Met():
    TLorentzVector(),
    MetUnclustEnUpDeltaX_(0),
    MetUnclustEnUpDeltaY_(0),
    covXX_(0),
    covXY_(0),
    covYY_(0),
    significance_(0),
    sumEt_(0),
    fiducialGenPhi_(0),
    fiducialGenPt_(0)
  {}
  Float_t MetUnclustEnUpDeltaX() const {return MetUnclustEnUpDeltaX_;}
  Float_t MetUnclustEnUpDeltaY() const {return MetUnclustEnUpDeltaY_;}
  Float_t covXX() const {return covXX_;}
  Float_t covXY() const {return covXY_;}
  Float_t covYY() const {return covYY_;}
  Float_t significance() const {return significance_;}
  Float_t sumEt() const {return sumEt_;}
  Float_t fiducialGenPhi() const {return fiducialGenPhi_;}
  Float_t fiducialGenPt() const {return fiducialGenPt_;}
private:
  Float_t MetUnclustEnUpDeltaX_;
  Float_t MetUnclustEnUpDeltaY_;
  Float_t covXX_;
  Float_t covXY_;
  Float_t covYY_;
  Float_t significance_;
  Float_t sumEt_;
  Float_t fiducialGenPhi_;
  Float_t fiducialGenPt_;
  void setMetUnclustEnUpDeltaX(const Float_t value) {MetUnclustEnUpDeltaX_ = value;}
  void setMetUnclustEnUpDeltaY(const Float_t value) {MetUnclustEnUpDeltaY_ = value;}
  void setcovXX(const Float_t value) {covXX_ = value;}
  void setcovXY(const Float_t value) {covXY_ = value;}
  void setcovYY(const Float_t value) {covYY_ = value;}
  void setsignificance(const Float_t value) {significance_ = value;}
  void setsumEt(const Float_t value) {sumEt_ = value;}
  void setfiducialGenPhi(const Float_t value) {fiducialGenPhi_ = value;}
  void setfiducialGenPt(const Float_t value) {fiducialGenPt_ = value;}
  void setLotentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Genmet: public TLorentzVector{
friend class URStreamer;
public:
//  Genmet():
//    
//  {}
  Genmet():
    TLorentzVector(),
    
  {}
  
private:
  
  
  void setLotentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Fatje: public TLorentzVector{
friend class URStreamer;
public:
//  Fatje(const Float_ &i_area_,const Float_ &i_btagCMVA_,const Float_ &i_btagCSVV2_,const Float_ &i_btagDeepB_,const Float_ &i_btagHbb_,const Float_ &i_msoftdrop_,const Float_ &i_msoftdrop_,const Float_ &i_n2b1_,const Float_ &i_n3b1_,const Float_ &i_tau1_,const Float_ &i_tau2_,const Float_ &i_tau3_,const Float_ &i_tau4_,const Int_ &i_jetId_,const Int_ &i_subJetIdx1_,const Int_ &i_subJetIdx2_):
//    
//  {}
  Fatje():
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
  Float_ area() const {return area_;}
  Float_ btagCMVA() const {return btagCMVA_;}
  Float_ btagCSVV2() const {return btagCSVV2_;}
  Float_ btagDeepB() const {return btagDeepB_;}
  Float_ btagHbb() const {return btagHbb_;}
  Float_ msoftdrop() const {return msoftdrop_;}
  Float_ msoftdrop() const {return msoftdrop_;}
  Float_ n2b1() const {return n2b1_;}
  Float_ n3b1() const {return n3b1_;}
  Float_ tau1() const {return tau1_;}
  Float_ tau2() const {return tau2_;}
  Float_ tau3() const {return tau3_;}
  Float_ tau4() const {return tau4_;}
  Int_ jetId() const {return jetId_;}
  Int_ subJetIdx1() const {return subJetIdx1_;}
  Int_ subJetIdx2() const {return subJetIdx2_;}
private:
  Float_ area_;
  Float_ btagCMVA_;
  Float_ btagCSVV2_;
  Float_ btagDeepB_;
  Float_ btagHbb_;
  Float_ msoftdrop_;
  Float_ msoftdrop_;
  Float_ n2b1_;
  Float_ n3b1_;
  Float_ tau1_;
  Float_ tau2_;
  Float_ tau3_;
  Float_ tau4_;
  Int_ jetId_;
  Int_ subJetIdx1_;
  Int_ subJetIdx2_;
  void setarea(const Float_ value) {area_ = value;}
  void setbtagCMVA(const Float_ value) {btagCMVA_ = value;}
  void setbtagCSVV2(const Float_ value) {btagCSVV2_ = value;}
  void setbtagDeepB(const Float_ value) {btagDeepB_ = value;}
  void setbtagHbb(const Float_ value) {btagHbb_ = value;}
  void setmsoftdrop(const Float_ value) {msoftdrop_ = value;}
  void setmsoftdrop(const Float_ value) {msoftdrop_ = value;}
  void setn2b1(const Float_ value) {n2b1_ = value;}
  void setn3b1(const Float_ value) {n3b1_ = value;}
  void settau1(const Float_ value) {tau1_ = value;}
  void settau2(const Float_ value) {tau2_ = value;}
  void settau3(const Float_ value) {tau3_ = value;}
  void settau4(const Float_ value) {tau4_ = value;}
  void setjetId(const Int_ value) {jetId_ = value;}
  void setsubJetIdx1(const Int_ value) {subJetIdx1_ = value;}
  void setsubJetIdx2(const Int_ value) {subJetIdx2_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Subje: public TLorentzVector{
friend class URStreamer;
public:
//  Subje(const Float_ &i_btagCMVA_,const Float_ &i_btagCSVV2_,const Float_ &i_btagDeepB_,const Float_ &i_n2b1_,const Float_ &i_n3b1_,const Float_ &i_tau1_,const Float_ &i_tau2_,const Float_ &i_tau3_,const Float_ &i_tau4_):
//    
//  {}
  Subje():
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
  Float_ btagCMVA() const {return btagCMVA_;}
  Float_ btagCSVV2() const {return btagCSVV2_;}
  Float_ btagDeepB() const {return btagDeepB_;}
  Float_ n2b1() const {return n2b1_;}
  Float_ n3b1() const {return n3b1_;}
  Float_ tau1() const {return tau1_;}
  Float_ tau2() const {return tau2_;}
  Float_ tau3() const {return tau3_;}
  Float_ tau4() const {return tau4_;}
private:
  Float_ btagCMVA_;
  Float_ btagCSVV2_;
  Float_ btagDeepB_;
  Float_ n2b1_;
  Float_ n3b1_;
  Float_ tau1_;
  Float_ tau2_;
  Float_ tau3_;
  Float_ tau4_;
  void setbtagCMVA(const Float_ value) {btagCMVA_ = value;}
  void setbtagCSVV2(const Float_ value) {btagCSVV2_ = value;}
  void setbtagDeepB(const Float_ value) {btagDeepB_ = value;}
  void setn2b1(const Float_ value) {n2b1_ = value;}
  void setn3b1(const Float_ value) {n3b1_ = value;}
  void settau1(const Float_ value) {tau1_ = value;}
  void settau2(const Float_ value) {tau2_ = value;}
  void settau3(const Float_ value) {tau3_ = value;}
  void settau4(const Float_ value) {tau4_ = value;}
  void setLotentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Flag{
friend class URStreamer;
public:
//  Flag(const Bool_t &i_HBHENoiseFilter_,const Bool_t &i_HBHENoiseIsoFilter_,const Bool_t &i_CSCTightHaloFilter_,const Bool_t &i_CSCTightHaloTrkMuUnvetoFilter_,const Bool_t &i_CSCTightHalo2015Filter_,const Bool_t &i_globalTightHalo2016Filter_,const Bool_t &i_globalSuperTightHalo2016Filter_,const Bool_t &i_HcalStripHaloFilter_,const Bool_t &i_hcalLaserEventFilter_,const Bool_t &i_EcalDeadCellTriggerPrimitiveFilter_,const Bool_t &i_EcalDeadCellBoundaryEnergyFilter_,const Bool_t &i_goodVertices_,const Bool_t &i_eeBadScFilter_,const Bool_t &i_ecalLaserCorrFilter_,const Bool_t &i_trkPOGFilters_,const Bool_t &i_chargedHadronTrackResolutionFilter_,const Bool_t &i_muonBadTrackFilter_,const Bool_t &i_trkPOG_,const Bool_t &i_trkPOG_,const Bool_t &i_trkPOG_,const Bool_t &i_METFilters_):
//    
//  {}
  Flag():
    HBHENoiseFilter_(0),
    HBHENoiseIsoFilter_(0),
    CSCTightHaloFilter_(0),
    CSCTightHaloTrkMuUnvetoFilter_(0),
    CSCTightHalo2015Filter_(0),
    globalTightHalo2016Filter_(0),
    globalSuperTightHalo2016Filter_(0),
    HcalStripHaloFilter_(0),
    hcalLaserEventFilter_(0),
    EcalDeadCellTriggerPrimitiveFilter_(0),
    EcalDeadCellBoundaryEnergyFilter_(0),
    goodVertices_(0),
    eeBadScFilter_(0),
    ecalLaserCorrFilter_(0),
    trkPOGFilters_(0),
    chargedHadronTrackResolutionFilter_(0),
    muonBadTrackFilter_(0),
    trkPOG_(0),
    trkPOG_(0),
    trkPOG_(0),
    METFilters_(0)
  {}
  Bool_t HBHENoiseFilter() const {return HBHENoiseFilter_;}
  Bool_t HBHENoiseIsoFilter() const {return HBHENoiseIsoFilter_;}
  Bool_t CSCTightHaloFilter() const {return CSCTightHaloFilter_;}
  Bool_t CSCTightHaloTrkMuUnvetoFilter() const {return CSCTightHaloTrkMuUnvetoFilter_;}
  Bool_t CSCTightHalo2015Filter() const {return CSCTightHalo2015Filter_;}
  Bool_t globalTightHalo2016Filter() const {return globalTightHalo2016Filter_;}
  Bool_t globalSuperTightHalo2016Filter() const {return globalSuperTightHalo2016Filter_;}
  Bool_t HcalStripHaloFilter() const {return HcalStripHaloFilter_;}
  Bool_t hcalLaserEventFilter() const {return hcalLaserEventFilter_;}
  Bool_t EcalDeadCellTriggerPrimitiveFilter() const {return EcalDeadCellTriggerPrimitiveFilter_;}
  Bool_t EcalDeadCellBoundaryEnergyFilter() const {return EcalDeadCellBoundaryEnergyFilter_;}
  Bool_t goodVertices() const {return goodVertices_;}
  Bool_t eeBadScFilter() const {return eeBadScFilter_;}
  Bool_t ecalLaserCorrFilter() const {return ecalLaserCorrFilter_;}
  Bool_t trkPOGFilters() const {return trkPOGFilters_;}
  Bool_t chargedHadronTrackResolutionFilter() const {return chargedHadronTrackResolutionFilter_;}
  Bool_t muonBadTrackFilter() const {return muonBadTrackFilter_;}
  Bool_t trkPOG() const {return trkPOG_;}
  Bool_t trkPOG() const {return trkPOG_;}
  Bool_t trkPOG() const {return trkPOG_;}
  Bool_t METFilters() const {return METFilters_;}
private:
  Bool_t HBHENoiseFilter_;
  Bool_t HBHENoiseIsoFilter_;
  Bool_t CSCTightHaloFilter_;
  Bool_t CSCTightHaloTrkMuUnvetoFilter_;
  Bool_t CSCTightHalo2015Filter_;
  Bool_t globalTightHalo2016Filter_;
  Bool_t globalSuperTightHalo2016Filter_;
  Bool_t HcalStripHaloFilter_;
  Bool_t hcalLaserEventFilter_;
  Bool_t EcalDeadCellTriggerPrimitiveFilter_;
  Bool_t EcalDeadCellBoundaryEnergyFilter_;
  Bool_t goodVertices_;
  Bool_t eeBadScFilter_;
  Bool_t ecalLaserCorrFilter_;
  Bool_t trkPOGFilters_;
  Bool_t chargedHadronTrackResolutionFilter_;
  Bool_t muonBadTrackFilter_;
  Bool_t trkPOG_;
  Bool_t trkPOG_;
  Bool_t trkPOG_;
  Bool_t METFilters_;
  void setHBHENoiseFilter(const Bool_t value) {HBHENoiseFilter_ = value;}
  void setHBHENoiseIsoFilter(const Bool_t value) {HBHENoiseIsoFilter_ = value;}
  void setCSCTightHaloFilter(const Bool_t value) {CSCTightHaloFilter_ = value;}
  void setCSCTightHaloTrkMuUnvetoFilter(const Bool_t value) {CSCTightHaloTrkMuUnvetoFilter_ = value;}
  void setCSCTightHalo2015Filter(const Bool_t value) {CSCTightHalo2015Filter_ = value;}
  void setglobalTightHalo2016Filter(const Bool_t value) {globalTightHalo2016Filter_ = value;}
  void setglobalSuperTightHalo2016Filter(const Bool_t value) {globalSuperTightHalo2016Filter_ = value;}
  void setHcalStripHaloFilter(const Bool_t value) {HcalStripHaloFilter_ = value;}
  void sethcalLaserEventFilter(const Bool_t value) {hcalLaserEventFilter_ = value;}
  void setEcalDeadCellTriggerPrimitiveFilter(const Bool_t value) {EcalDeadCellTriggerPrimitiveFilter_ = value;}
  void setEcalDeadCellBoundaryEnergyFilter(const Bool_t value) {EcalDeadCellBoundaryEnergyFilter_ = value;}
  void setgoodVertices(const Bool_t value) {goodVertices_ = value;}
  void seteeBadScFilter(const Bool_t value) {eeBadScFilter_ = value;}
  void setecalLaserCorrFilter(const Bool_t value) {ecalLaserCorrFilter_ = value;}
  void settrkPOGFilters(const Bool_t value) {trkPOGFilters_ = value;}
  void setchargedHadronTrackResolutionFilter(const Bool_t value) {chargedHadronTrackResolutionFilter_ = value;}
  void setmuonBadTrackFilter(const Bool_t value) {muonBadTrackFilter_ = value;}
  void settrkPOG(const Bool_t value) {trkPOG_ = value;}
  void settrkPOG(const Bool_t value) {trkPOG_ = value;}
  void settrkPOG(const Bool_t value) {trkPOG_ = value;}
  void setMETFilters(const Bool_t value) {METFilters_ = value;}
};

class Pileup{
friend class URStreamer;
public:
//  Pileup(const Float_t &i_nTrueInt_,const Int_t &i_nPU_,const Int_t &i_sumEOOT_,const Int_t &i_sumLOOT_):
//    
//  {}
  Pileup():
    nTrueInt_(0),
    nPU_(0),
    sumEOOT_(0),
    sumLOOT_(0)
  {}
  Float_t nTrueInt() const {return nTrueInt_;}
  Int_t nPU() const {return nPU_;}
  Int_t sumEOOT() const {return sumEOOT_;}
  Int_t sumLOOT() const {return sumLOOT_;}
private:
  Float_t nTrueInt_;
  Int_t nPU_;
  Int_t sumEOOT_;
  Int_t sumLOOT_;
  void setnTrueInt(const Float_t value) {nTrueInt_ = value;}
  void setnPU(const Int_t value) {nPU_ = value;}
  void setsumEOOT(const Int_t value) {sumEOOT_ = value;}
  void setsumLOOT(const Int_t value) {sumLOOT_ = value;}
};

class Lheweight{
friend class URStreamer;
public:
//  Lheweight(const Float_t &i_originalXWGTUP_):
//    
//  {}
  Lheweight():
    originalXWGTUP_(0)
  {}
  Float_t originalXWGTUP() const {return originalXWGTUP_;}
private:
  Float_t originalXWGTUP_;
  void setoriginalXWGTUP(const Float_t value) {originalXWGTUP_ = value;}
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
  Float_t *LHEPdfWeight;
  UInt_t nLHEScaleWeight;
  Float_t *LHEScaleWeight;
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
    are_Jet_loaded_(0), Jet_(),
    are_GenJetAK8_loaded_(0), GenJetAK8_(),
    are_GenVisTau_loaded_(0), GenVisTau_(),
    are_CaloMET_loaded_(0), CaloMET_(),
    are_GenDressedLepton_loaded_(0), GenDressedLepton_(),
    are_PV_loaded_(0), PV_(),
    are_Generator_loaded_(0), Generator_(),
    are_TrigObj_loaded_(0), TrigObj_(),
    are_Photon_loaded_(0), Photon_(),
    are_GenJet_loaded_(0), GenJet_(),
    are_RawMET_loaded_(0), RawMET_(),
    are_Electron_loaded_(0), Electron_(),
    are_SoftActivityJet_loaded_(0), SoftActivityJet_(),
    are_L1simulation_loaded_(0), L1simulation_(),
    are_GenPart_loaded_(0), GenPart_(),
    are_LHE_loaded_(0), LHE_(),
    are_TkMET_loaded_(0), TkMET_(),
    are_Tau_loaded_(0), Tau_(),
    are_PuppiMET_loaded_(0), PuppiMET_(),
    are_Muon_loaded_(0), Muon_(),
    are_OtherPV_loaded_(0), OtherPV_(),
    are_HLT_loaded_(0), HLT_(),
    are_SV_loaded_(0), SV_(),
    are_MET_loaded_(0), MET_(),
    are_GenMET_loaded_(0), GenMET_(),
    are_FatJet_loaded_(0), FatJet_(),
    are_SubJet_loaded_(0), SubJet_(),
    are_Flag_loaded_(0), Flag_(),
    are_Pileup_loaded_(0), Pileup_(),
    are_LHEWeight_loaded_(0), LHEWeight_()
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
    Jet_.clear();
    GenJetAK8_.clear();
    GenVisTau_.clear();
    
    GenDressedLepton_.clear();
    
    
    TrigObj_.clear();
    Photon_.clear();
    GenJet_.clear();
    
    Electron_.clear();
    
    
    GenPart_.clear();
    
    
    Tau_.clear();
    
    Muon_.clear();
    OtherPV_.clear();
    
    SV_.clear();
    
    
    FatJet_.clear();
    SubJet_.clear();
    
    
    
    current_entry_++;
    if(current_entry_ < entries_){
      tree_->GetEntry(current_entry_);
      return true;
    }
    return false;
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
  
  void loadCalomet(){
    if(!are_CaloMET_loaded_){
      tree_->SetBranchStatus("CaloMET_phi", 1); tree_->SetBranchAddress("CaloMET_phi", &CaloMET_phi_);
      tree_->SetBranchStatus("CaloMET_pt", 1); tree_->SetBranchAddress("CaloMET_pt", &CaloMET_pt_);
      tree_->SetBranchStatus("CaloMET_sumEt", 1); tree_->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt_);
      are_CaloMET_loaded_ = true;
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
  
  void loadPv(){
    if(!are_PV_loaded_){
      tree_->SetBranchStatus("PV_ndof", 1); tree_->SetBranchAddress("PV_ndof", &PV_ndof_);
      tree_->SetBranchStatus("PV_x", 1); tree_->SetBranchAddress("PV_x", &PV_x_);
      tree_->SetBranchStatus("PV_y", 1); tree_->SetBranchAddress("PV_y", &PV_y_);
      tree_->SetBranchStatus("PV_z", 1); tree_->SetBranchAddress("PV_z", &PV_z_);
      tree_->SetBranchStatus("PV_chi2", 1); tree_->SetBranchAddress("PV_chi2", &PV_chi2_);
      tree_->SetBranchStatus("PV_score", 1); tree_->SetBranchAddress("PV_score", &PV_score_);
      tree_->SetBranchStatus("PV_npvs", 1); tree_->SetBranchAddress("PV_npvs", &PV_npvs_);
      tree_->SetBranchStatus("PV_npvsGood", 1); tree_->SetBranchAddress("PV_npvsGood", &PV_npvsGood_);
      are_PV_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenerator(){
    if(!are_Generator_loaded_){
      tree_->SetBranchStatus("Generator_binvar", 1); tree_->SetBranchAddress("Generator_binvar", &Generator_binvar_);
      tree_->SetBranchStatus("Generator_scalePDF", 1); tree_->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF_);
      tree_->SetBranchStatus("Generator_weight", 1); tree_->SetBranchAddress("Generator_weight", &Generator_weight_);
      tree_->SetBranchStatus("Generator_x1", 1); tree_->SetBranchAddress("Generator_x1", &Generator_x1_);
      tree_->SetBranchStatus("Generator_x2", 1); tree_->SetBranchAddress("Generator_x2", &Generator_x2_);
      tree_->SetBranchStatus("Generator_xpdf1", 1); tree_->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1_);
      tree_->SetBranchStatus("Generator_xpdf2", 1); tree_->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2_);
      tree_->SetBranchStatus("Generator_id1", 1); tree_->SetBranchAddress("Generator_id1", &Generator_id1_);
      tree_->SetBranchStatus("Generator_id2", 1); tree_->SetBranchAddress("Generator_id2", &Generator_id2_);
      are_Generator_loaded_ = true;
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
  
  void loadRawmet(){
    if(!are_RawMET_loaded_){
      tree_->SetBranchStatus("RawMET_phi", 1); tree_->SetBranchAddress("RawMET_phi", &RawMET_phi_);
      tree_->SetBranchStatus("RawMET_pt", 1); tree_->SetBranchAddress("RawMET_pt", &RawMET_pt_);
      tree_->SetBranchStatus("RawMET_sumEt", 1); tree_->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt_);
      are_RawMET_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
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
  
  void loadL1Simulation(){
    if(!are_L1simulation_loaded_){
      tree_->SetBranchStatus("L1simulation_step", 1); tree_->SetBranchAddress("L1simulation_step", &L1simulation_step_);
      are_L1simulation_loaded_ = true;
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
  
  void loadLhe(){
    if(!are_LHE_loaded_){
      tree_->SetBranchStatus("LHEWeight_originalXWGTUP", 1); tree_->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP_);
      tree_->SetBranchStatus("LHEPdfWeight", 1); tree_->SetBranchAddress("LHEPdfWeight", &LHEPdfWeight);
      tree_->SetBranchStatus("LHEScaleWeight", 1); tree_->SetBranchAddress("LHEScaleWeight", &LHEScaleWeight);
      tree_->SetBranchStatus("LHE_HT", 1); tree_->SetBranchAddress("LHE_HT", &LHE_HT_);
      tree_->SetBranchStatus("LHE_HTIncoming", 1); tree_->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming_);
      tree_->SetBranchStatus("LHE_Vpt", 1); tree_->SetBranchAddress("LHE_Vpt", &LHE_Vpt_);
      tree_->SetBranchStatus("LHE_Njets", 1); tree_->SetBranchAddress("LHE_Njets", &LHE_Njets_);
      tree_->SetBranchStatus("LHE_Nb", 1); tree_->SetBranchAddress("LHE_Nb", &LHE_Nb_);
      tree_->SetBranchStatus("LHE_Nc", 1); tree_->SetBranchAddress("LHE_Nc", &LHE_Nc_);
      tree_->SetBranchStatus("LHE_Nuds", 1); tree_->SetBranchAddress("LHE_Nuds", &LHE_Nuds_);
      tree_->SetBranchStatus("LHE_Nglu", 1); tree_->SetBranchAddress("LHE_Nglu", &LHE_Nglu_);
      tree_->SetBranchStatus("LHE_NpNLO", 1); tree_->SetBranchAddress("LHE_NpNLO", &LHE_NpNLO_);
      tree_->SetBranchStatus("LHE_NpLO", 1); tree_->SetBranchAddress("LHE_NpLO", &LHE_NpLO_);
      are_LHE_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadTkmet(){
    if(!are_TkMET_loaded_){
      tree_->SetBranchStatus("TkMET_phi", 1); tree_->SetBranchAddress("TkMET_phi", &TkMET_phi_);
      tree_->SetBranchStatus("TkMET_pt", 1); tree_->SetBranchAddress("TkMET_pt", &TkMET_pt_);
      tree_->SetBranchStatus("TkMET_sumEt", 1); tree_->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt_);
      are_TkMET_loaded_ = true;
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
  
  void loadPuppimet(){
    if(!are_PuppiMET_loaded_){
      tree_->SetBranchStatus("PuppiMET_phi", 1); tree_->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi_);
      tree_->SetBranchStatus("PuppiMET_pt", 1); tree_->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt_);
      tree_->SetBranchStatus("PuppiMET_sumEt", 1); tree_->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt_);
      are_PuppiMET_loaded_ = true;
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
  
  void loadOtherpv(){
    if(!are_OtherPV_loaded_){
      tree_->SetBranchStatus("OtherPV_z", 1); tree_->SetBranchAddress("OtherPV_z", &OtherPV_z_);
      are_OtherPV_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadHlt(){
    if(!are_HLT_loaded_){
      tree_->SetBranchStatus("HLTriggerFirstPath", 1); tree_->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath);
      tree_->SetBranchStatus("HLT_AK8PFJet360_TrimMass30", 1); tree_->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30_);
      tree_->SetBranchStatus("HLT_AK8PFJet400_TrimMass30", 1); tree_->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30_);
      tree_->SetBranchStatus("HLT_AK8PFHT750_TrimMass50", 1); tree_->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50_);
      tree_->SetBranchStatus("HLT_AK8PFHT800_TrimMass50", 1); tree_->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50_);
      tree_->SetBranchStatus("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20", 1); tree_->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_);
      tree_->SetBranchStatus("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087", 1); tree_->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_);
      tree_->SetBranchStatus("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087", 1); tree_->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_);
      tree_->SetBranchStatus("HLT_AK8DiPFJet300_200_TrimMass30", 1); tree_->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30", &HLT_AK8DiPFJet300_200_TrimMass30_);
      tree_->SetBranchStatus("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", 1); tree_->SetBranchAddress("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_);
      tree_->SetBranchStatus("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50", 1); tree_->SetBranchAddress("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_);
      tree_->SetBranchStatus("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20", 1); tree_->SetBranchAddress("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20", &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_);
      tree_->SetBranchStatus("HLT_AK8DiPFJet280_200_TrimMass30", 1); tree_->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30", &HLT_AK8DiPFJet280_200_TrimMass30_);
      tree_->SetBranchStatus("HLT_AK8DiPFJet250_200_TrimMass30", 1); tree_->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30", &HLT_AK8DiPFJet250_200_TrimMass30_);
      tree_->SetBranchStatus("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20", 1); tree_->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_);
      tree_->SetBranchStatus("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20", 1); tree_->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_);
      tree_->SetBranchStatus("HLT_CaloJet260", 1); tree_->SetBranchAddress("HLT_CaloJet260", &HLT_CaloJet260_);
      tree_->SetBranchStatus("HLT_CaloJet500_NoJetID", 1); tree_->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID_);
      tree_->SetBranchStatus("HLT_Dimuon13_PsiPrime", 1); tree_->SetBranchAddress("HLT_Dimuon13_PsiPrime", &HLT_Dimuon13_PsiPrime_);
      tree_->SetBranchStatus("HLT_Dimuon13_Upsilon", 1); tree_->SetBranchAddress("HLT_Dimuon13_Upsilon", &HLT_Dimuon13_Upsilon_);
      tree_->SetBranchStatus("HLT_Dimuon20_Jpsi", 1); tree_->SetBranchAddress("HLT_Dimuon20_Jpsi", &HLT_Dimuon20_Jpsi_);
      tree_->SetBranchStatus("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf", &HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_DoubleEle25_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("HLT_DoubleEle25_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle25_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("HLT_DoubleEle33_CaloIdL", 1); tree_->SetBranchAddress("HLT_DoubleEle33_CaloIdL", &HLT_DoubleEle33_CaloIdL_);
      tree_->SetBranchStatus("HLT_DoubleEle33_CaloIdL_MW", 1); tree_->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW_);
      tree_->SetBranchStatus("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW", 1); tree_->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_);
      tree_->SetBranchStatus("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1", 1); tree_->SetBranchAddress("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1", &HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_);
      tree_->SetBranchStatus("HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1", 1); tree_->SetBranchAddress("HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1", &HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_);
      tree_->SetBranchStatus("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1", 1); tree_->SetBranchAddress("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1", &HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_);
      tree_->SetBranchStatus("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("HLT_DoubleMu33NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_DoubleMu33NoFiltersNoVtx", &HLT_DoubleMu33NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_DoubleMu38NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_DoubleMu38NoFiltersNoVtx", &HLT_DoubleMu38NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_DoubleMu23NoFiltersNoVtxDisplaced", 1); tree_->SetBranchAddress("HLT_DoubleMu23NoFiltersNoVtxDisplaced", &HLT_DoubleMu23NoFiltersNoVtxDisplaced_);
      tree_->SetBranchStatus("HLT_DoubleMu28NoFiltersNoVtxDisplaced", 1); tree_->SetBranchAddress("HLT_DoubleMu28NoFiltersNoVtxDisplaced", &HLT_DoubleMu28NoFiltersNoVtxDisplaced_);
      tree_->SetBranchStatus("HLT_DoubleMu0", 1); tree_->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0_);
      tree_->SetBranchStatus("HLT_DoubleMu4_3_Bs", 1); tree_->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs_);
      tree_->SetBranchStatus("HLT_DoubleMu4_3_Jpsi_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_Displaced", &HLT_DoubleMu4_3_Jpsi_Displaced_);
      tree_->SetBranchStatus("HLT_DoubleMu4_JpsiTrk_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced_);
      tree_->SetBranchStatus("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_);
      tree_->SetBranchStatus("HLT_DoubleMu3_Trk_Tau3mu", 1); tree_->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu_);
      tree_->SetBranchStatus("HLT_DoubleMu4_PsiPrimeTrk_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced_);
      tree_->SetBranchStatus("HLT_Mu7p5_L2Mu2_Jpsi", 1); tree_->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi_);
      tree_->SetBranchStatus("HLT_Mu7p5_L2Mu2_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track2_Jpsi", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track3p5_Jpsi", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track7_Jpsi", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track2_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track3p5_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track7_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon_);
      tree_->SetBranchStatus("HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing", 1); tree_->SetBranchAddress("HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing", &HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_);
      tree_->SetBranchStatus("HLT_Dimuon0er16_Jpsi_NoVertexing", 1); tree_->SetBranchAddress("HLT_Dimuon0er16_Jpsi_NoVertexing", &HLT_Dimuon0er16_Jpsi_NoVertexing_);
      tree_->SetBranchStatus("HLT_Dimuon6_Jpsi_NoVertexing", 1); tree_->SetBranchAddress("HLT_Dimuon6_Jpsi_NoVertexing", &HLT_Dimuon6_Jpsi_NoVertexing_);
      tree_->SetBranchStatus("HLT_Photon150", 1); tree_->SetBranchAddress("HLT_Photon150", &HLT_Photon150_);
      tree_->SetBranchStatus("HLT_Photon90_CaloIdL_HT300", 1); tree_->SetBranchAddress("HLT_Photon90_CaloIdL_HT300", &HLT_Photon90_CaloIdL_HT300_);
      tree_->SetBranchStatus("HLT_HT250_CaloMET70", 1); tree_->SetBranchAddress("HLT_HT250_CaloMET70", &HLT_HT250_CaloMET70_);
      tree_->SetBranchStatus("HLT_DoublePhoton60", 1); tree_->SetBranchAddress("HLT_DoublePhoton60", &HLT_DoublePhoton60_);
      tree_->SetBranchStatus("HLT_DoublePhoton85", 1); tree_->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85_);
      tree_->SetBranchStatus("HLT_Ele17_Ele8_Gsf", 1); tree_->SetBranchAddress("HLT_Ele17_Ele8_Gsf", &HLT_Ele17_Ele8_Gsf_);
      tree_->SetBranchStatus("HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28", 1); tree_->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28", &HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28_);
      tree_->SetBranchStatus("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29", 1); tree_->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29", &HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29_);
      tree_->SetBranchStatus("HLT_Ele22_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf", &HLT_Ele22_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", 1); tree_->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
      tree_->SetBranchStatus("HLT_Ele23_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele23_WPLoose_Gsf", &HLT_Ele23_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele23_WPLoose_Gsf_WHbbBoost", 1); tree_->SetBranchAddress("HLT_Ele23_WPLoose_Gsf_WHbbBoost", &HLT_Ele23_WPLoose_Gsf_WHbbBoost_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf", &HLT_Ele24_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_);
      tree_->SetBranchStatus("HLT_Ele25_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele25_WPTight_Gsf", &HLT_Ele25_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele25_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele25_eta2p1_WPLoose_Gsf", &HLT_Ele25_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele25_eta2p1_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf", &HLT_Ele25_eta2p1_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele27_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele27_WPLoose_Gsf", &HLT_Ele27_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele27_WPLoose_Gsf_WHbbBoost", 1); tree_->SetBranchAddress("HLT_Ele27_WPLoose_Gsf_WHbbBoost", &HLT_Ele27_WPLoose_Gsf_WHbbBoost_);
      tree_->SetBranchStatus("HLT_Ele27_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele27_WPTight_Gsf_L1JetTauSeeded", 1); tree_->SetBranchAddress("HLT_Ele27_WPTight_Gsf_L1JetTauSeeded", &HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_);
      tree_->SetBranchStatus("HLT_Ele27_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf", &HLT_Ele27_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", 1); tree_->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
      tree_->SetBranchStatus("HLT_Ele27_eta2p1_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele27_eta2p1_WPTight_Gsf", &HLT_Ele27_eta2p1_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele30_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele30_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele30_eta2p1_WPLoose_Gsf", &HLT_Ele30_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele30_eta2p1_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf", &HLT_Ele30_eta2p1_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele32_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele32_eta2p1_WPLoose_Gsf", &HLT_Ele32_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", 1); tree_->SetBranchAddress("HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
      tree_->SetBranchStatus("HLT_Ele32_eta2p1_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf", &HLT_Ele32_eta2p1_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele35_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele35_WPLoose_Gsf", &HLT_Ele35_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50", 1); tree_->SetBranchAddress("HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50", &HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_);
      tree_->SetBranchStatus("HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", 1); tree_->SetBranchAddress("HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
      tree_->SetBranchStatus("HLT_Ele45_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele45_WPLoose_Gsf", &HLT_Ele45_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded", 1); tree_->SetBranchAddress("HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded", &HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded_);
      tree_->SetBranchStatus("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", 1); tree_->SetBranchAddress("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", &HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_);
      tree_->SetBranchStatus("HLT_Ele105_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele105_CaloIdVT_GsfTrkIdT", &HLT_Ele105_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_Ele30WP60_SC4_Mass55", 1); tree_->SetBranchAddress("HLT_Ele30WP60_SC4_Mass55", &HLT_Ele30WP60_SC4_Mass55_);
      tree_->SetBranchStatus("HLT_Ele30WP60_Ele8_Mass55", 1); tree_->SetBranchAddress("HLT_Ele30WP60_Ele8_Mass55", &HLT_Ele30WP60_Ele8_Mass55_);
      tree_->SetBranchStatus("HLT_HT200", 1); tree_->SetBranchAddress("HLT_HT200", &HLT_HT200_);
      tree_->SetBranchStatus("HLT_HT275", 1); tree_->SetBranchAddress("HLT_HT275", &HLT_HT275_);
      tree_->SetBranchStatus("HLT_HT325", 1); tree_->SetBranchAddress("HLT_HT325", &HLT_HT325_);
      tree_->SetBranchStatus("HLT_HT425", 1); tree_->SetBranchAddress("HLT_HT425", &HLT_HT425_);
      tree_->SetBranchStatus("HLT_HT575", 1); tree_->SetBranchAddress("HLT_HT575", &HLT_HT575_);
      tree_->SetBranchStatus("HLT_HT410to430", 1); tree_->SetBranchAddress("HLT_HT410to430", &HLT_HT410to430_);
      tree_->SetBranchStatus("HLT_HT430to450", 1); tree_->SetBranchAddress("HLT_HT430to450", &HLT_HT430to450_);
      tree_->SetBranchStatus("HLT_HT450to470", 1); tree_->SetBranchAddress("HLT_HT450to470", &HLT_HT450to470_);
      tree_->SetBranchStatus("HLT_HT470to500", 1); tree_->SetBranchAddress("HLT_HT470to500", &HLT_HT470to500_);
      tree_->SetBranchStatus("HLT_HT500to550", 1); tree_->SetBranchAddress("HLT_HT500to550", &HLT_HT500to550_);
      tree_->SetBranchStatus("HLT_HT550to650", 1); tree_->SetBranchAddress("HLT_HT550to650", &HLT_HT550to650_);
      tree_->SetBranchStatus("HLT_HT650", 1); tree_->SetBranchAddress("HLT_HT650", &HLT_HT650_);
      tree_->SetBranchStatus("HLT_Mu16_eta2p1_MET30", 1); tree_->SetBranchAddress("HLT_Mu16_eta2p1_MET30", &HLT_Mu16_eta2p1_MET30_);
      tree_->SetBranchStatus("HLT_IsoMu16_eta2p1_MET30", 1); tree_->SetBranchAddress("HLT_IsoMu16_eta2p1_MET30", &HLT_IsoMu16_eta2p1_MET30_);
      tree_->SetBranchStatus("HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1", 1); tree_->SetBranchAddress("HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1", &HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1_);
      tree_->SetBranchStatus("HLT_IsoMu17_eta2p1", 1); tree_->SetBranchAddress("HLT_IsoMu17_eta2p1", &HLT_IsoMu17_eta2p1_);
      tree_->SetBranchStatus("HLT_IsoMu17_eta2p1_LooseIsoPFTau20", 1); tree_->SetBranchAddress("HLT_IsoMu17_eta2p1_LooseIsoPFTau20", &HLT_IsoMu17_eta2p1_LooseIsoPFTau20_);
      tree_->SetBranchStatus("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1", 1); tree_->SetBranchAddress("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_);
      tree_->SetBranchStatus("HLT_DoubleIsoMu17_eta2p1", 1); tree_->SetBranchAddress("HLT_DoubleIsoMu17_eta2p1", &HLT_DoubleIsoMu17_eta2p1_);
      tree_->SetBranchStatus("HLT_DoubleIsoMu17_eta2p1_noDzCut", 1); tree_->SetBranchAddress("HLT_DoubleIsoMu17_eta2p1_noDzCut", &HLT_DoubleIsoMu17_eta2p1_noDzCut_);
      tree_->SetBranchStatus("HLT_IsoMu18", 1); tree_->SetBranchAddress("HLT_IsoMu18", &HLT_IsoMu18_);
      tree_->SetBranchStatus("HLT_IsoMu19_eta2p1_LooseIsoPFTau20", 1); tree_->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseIsoPFTau20", &HLT_IsoMu19_eta2p1_LooseIsoPFTau20_);
      tree_->SetBranchStatus("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1", 1); tree_->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_);
      tree_->SetBranchStatus("HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20", 1); tree_->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20", &HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_);
      tree_->SetBranchStatus("HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_IsoMu20", 1); tree_->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20_);
      tree_->SetBranchStatus("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1", 1); tree_->SetBranchAddress("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_);
      tree_->SetBranchStatus("HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1", 1); tree_->SetBranchAddress("HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1", &HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1_);
      tree_->SetBranchStatus("HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_IsoMu22", 1); tree_->SetBranchAddress("HLT_IsoMu22", &HLT_IsoMu22_);
      tree_->SetBranchStatus("HLT_IsoMu22_eta2p1", 1); tree_->SetBranchAddress("HLT_IsoMu22_eta2p1", &HLT_IsoMu22_eta2p1_);
      tree_->SetBranchStatus("HLT_IsoMu24", 1); tree_->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24_);
      tree_->SetBranchStatus("HLT_IsoMu27", 1); tree_->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27_);
      tree_->SetBranchStatus("HLT_IsoTkMu18", 1); tree_->SetBranchAddress("HLT_IsoTkMu18", &HLT_IsoTkMu18_);
      tree_->SetBranchStatus("HLT_IsoTkMu20", 1); tree_->SetBranchAddress("HLT_IsoTkMu20", &HLT_IsoTkMu20_);
      tree_->SetBranchStatus("HLT_IsoTkMu22", 1); tree_->SetBranchAddress("HLT_IsoTkMu22", &HLT_IsoTkMu22_);
      tree_->SetBranchStatus("HLT_IsoTkMu22_eta2p1", 1); tree_->SetBranchAddress("HLT_IsoTkMu22_eta2p1", &HLT_IsoTkMu22_eta2p1_);
      tree_->SetBranchStatus("HLT_IsoTkMu24", 1); tree_->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24_);
      tree_->SetBranchStatus("HLT_IsoTkMu27", 1); tree_->SetBranchAddress("HLT_IsoTkMu27", &HLT_IsoTkMu27_);
      tree_->SetBranchStatus("HLT_JetE30_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_JetE30_NoBPTX3BX", &HLT_JetE30_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_JetE30_NoBPTX", 1); tree_->SetBranchAddress("HLT_JetE30_NoBPTX", &HLT_JetE30_NoBPTX_);
      tree_->SetBranchStatus("HLT_JetE50_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_JetE50_NoBPTX3BX", &HLT_JetE50_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_JetE70_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_JetE70_NoBPTX3BX", &HLT_JetE70_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_L1SingleMu18", 1); tree_->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18_);
      tree_->SetBranchStatus("HLT_L2Mu10", 1); tree_->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10_);
      tree_->SetBranchStatus("HLT_L1SingleMuOpen", 1); tree_->SetBranchAddress("HLT_L1SingleMuOpen", &HLT_L1SingleMuOpen_);
      tree_->SetBranchStatus("HLT_L1SingleMuOpen_DT", 1); tree_->SetBranchAddress("HLT_L1SingleMuOpen_DT", &HLT_L1SingleMuOpen_DT_);
      tree_->SetBranchStatus("HLT_L2DoubleMu23_NoVertex", 1); tree_->SetBranchAddress("HLT_L2DoubleMu23_NoVertex", &HLT_L2DoubleMu23_NoVertex_);
      tree_->SetBranchStatus("HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10", 1); tree_->SetBranchAddress("HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10", &HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_);
      tree_->SetBranchStatus("HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10", 1); tree_->SetBranchAddress("HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10", &HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10_);
      tree_->SetBranchStatus("HLT_L2Mu10_NoVertex_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_L2Mu10_NoVertex_NoBPTX", 1); tree_->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX_);
      tree_->SetBranchStatus("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_LooseIsoPFTau50_Trk30_eta2p1", 1); tree_->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1", &HLT_LooseIsoPFTau50_Trk30_eta2p1_);
      tree_->SetBranchStatus("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80", 1); tree_->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_);
      tree_->SetBranchStatus("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90", 1); tree_->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_);
      tree_->SetBranchStatus("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110", 1); tree_->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_);
      tree_->SetBranchStatus("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120", 1); tree_->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_);
      tree_->SetBranchStatus("HLT_PFTau120_eta2p1", 1); tree_->SetBranchAddress("HLT_PFTau120_eta2p1", &HLT_PFTau120_eta2p1_);
      tree_->SetBranchStatus("HLT_PFTau140_eta2p1", 1); tree_->SetBranchAddress("HLT_PFTau140_eta2p1", &HLT_PFTau140_eta2p1_);
      tree_->SetBranchStatus("HLT_VLooseIsoPFTau120_Trk50_eta2p1", 1); tree_->SetBranchAddress("HLT_VLooseIsoPFTau120_Trk50_eta2p1", &HLT_VLooseIsoPFTau120_Trk50_eta2p1_);
      tree_->SetBranchStatus("HLT_VLooseIsoPFTau140_Trk50_eta2p1", 1); tree_->SetBranchAddress("HLT_VLooseIsoPFTau140_Trk50_eta2p1", &HLT_VLooseIsoPFTau140_Trk50_eta2p1_);
      tree_->SetBranchStatus("HLT_Mu17_Mu8", 1); tree_->SetBranchAddress("HLT_Mu17_Mu8", &HLT_Mu17_Mu8_);
      tree_->SetBranchStatus("HLT_Mu17_Mu8_DZ", 1); tree_->SetBranchAddress("HLT_Mu17_Mu8_DZ", &HLT_Mu17_Mu8_DZ_);
      tree_->SetBranchStatus("HLT_Mu17_Mu8_SameSign", 1); tree_->SetBranchAddress("HLT_Mu17_Mu8_SameSign", &HLT_Mu17_Mu8_SameSign_);
      tree_->SetBranchStatus("HLT_Mu17_Mu8_SameSign_DZ", 1); tree_->SetBranchAddress("HLT_Mu17_Mu8_SameSign_DZ", &HLT_Mu17_Mu8_SameSign_DZ_);
      tree_->SetBranchStatus("HLT_Mu20_Mu10", 1); tree_->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10_);
      tree_->SetBranchStatus("HLT_Mu20_Mu10_DZ", 1); tree_->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ_);
      tree_->SetBranchStatus("HLT_Mu20_Mu10_SameSign", 1); tree_->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign_);
      tree_->SetBranchStatus("HLT_Mu20_Mu10_SameSign_DZ", 1); tree_->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ_);
      tree_->SetBranchStatus("HLT_Mu17_TkMu8_DZ", 1); tree_->SetBranchAddress("HLT_Mu17_TkMu8_DZ", &HLT_Mu17_TkMu8_DZ_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu25_TkMu0_dEta18_Onia", 1); tree_->SetBranchAddress("HLT_Mu25_TkMu0_dEta18_Onia", &HLT_Mu25_TkMu0_dEta18_Onia_);
      tree_->SetBranchStatus("HLT_Mu27_TkMu8", 1); tree_->SetBranchAddress("HLT_Mu27_TkMu8", &HLT_Mu27_TkMu8_);
      tree_->SetBranchStatus("HLT_Mu30_TkMu11", 1); tree_->SetBranchAddress("HLT_Mu30_TkMu11", &HLT_Mu30_TkMu11_);
      tree_->SetBranchStatus("HLT_Mu30_eta2p1_PFJet150_PFJet50", 1); tree_->SetBranchAddress("HLT_Mu30_eta2p1_PFJet150_PFJet50", &HLT_Mu30_eta2p1_PFJet150_PFJet50_);
      tree_->SetBranchStatus("HLT_Mu40_TkMu11", 1); tree_->SetBranchAddress("HLT_Mu40_TkMu11", &HLT_Mu40_TkMu11_);
      tree_->SetBranchStatus("HLT_Mu40_eta2p1_PFJet200_PFJet50", 1); tree_->SetBranchAddress("HLT_Mu40_eta2p1_PFJet200_PFJet50", &HLT_Mu40_eta2p1_PFJet200_PFJet50_);
      tree_->SetBranchStatus("HLT_Mu20", 1); tree_->SetBranchAddress("HLT_Mu20", &HLT_Mu20_);
      tree_->SetBranchStatus("HLT_TkMu17", 1); tree_->SetBranchAddress("HLT_TkMu17", &HLT_TkMu17_);
      tree_->SetBranchStatus("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", 1); tree_->SetBranchAddress("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_);
      tree_->SetBranchStatus("HLT_TkMu20", 1); tree_->SetBranchAddress("HLT_TkMu20", &HLT_TkMu20_);
      tree_->SetBranchStatus("HLT_Mu24_eta2p1", 1); tree_->SetBranchAddress("HLT_Mu24_eta2p1", &HLT_Mu24_eta2p1_);
      tree_->SetBranchStatus("HLT_TkMu24_eta2p1", 1); tree_->SetBranchAddress("HLT_TkMu24_eta2p1", &HLT_TkMu24_eta2p1_);
      tree_->SetBranchStatus("HLT_Mu27", 1); tree_->SetBranchAddress("HLT_Mu27", &HLT_Mu27_);
      tree_->SetBranchStatus("HLT_TkMu27", 1); tree_->SetBranchAddress("HLT_TkMu27", &HLT_TkMu27_);
      tree_->SetBranchStatus("HLT_Mu45_eta2p1", 1); tree_->SetBranchAddress("HLT_Mu45_eta2p1", &HLT_Mu45_eta2p1_);
      tree_->SetBranchStatus("HLT_Mu50", 1); tree_->SetBranchAddress("HLT_Mu50", &HLT_Mu50_);
      tree_->SetBranchStatus("HLT_TkMu50", 1); tree_->SetBranchAddress("HLT_TkMu50", &HLT_TkMu50_);
      tree_->SetBranchStatus("HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL", &HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_);
      tree_->SetBranchStatus("HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL", &HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL_);
      tree_->SetBranchStatus("HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL", &HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL_);
      tree_->SetBranchStatus("HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL", &HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL_);
      tree_->SetBranchStatus("HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL", &HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL_);
      tree_->SetBranchStatus("HLT_DoubleMu18NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_DoubleMu18NoFiltersNoVtx", &HLT_DoubleMu18NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight", 1); tree_->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight", &HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight_);
      tree_->SetBranchStatus("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose", 1); tree_->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose", &HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose_);
      tree_->SetBranchStatus("HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose", 1); tree_->SetBranchAddress("HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose", &HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose_);
      tree_->SetBranchStatus("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight", 1); tree_->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight", &HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight_);
      tree_->SetBranchStatus("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose", 1); tree_->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose", &HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose_);
      tree_->SetBranchStatus("HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose", 1); tree_->SetBranchAddress("HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose", &HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose_);
      tree_->SetBranchStatus("HLT_Mu28NoFiltersNoVtx_CentralCaloJet40", 1); tree_->SetBranchAddress("HLT_Mu28NoFiltersNoVtx_CentralCaloJet40", &HLT_Mu28NoFiltersNoVtx_CentralCaloJet40_);
      tree_->SetBranchStatus("HLT_PFHT300_PFMET100", 1); tree_->SetBranchAddress("HLT_PFHT300_PFMET100", &HLT_PFHT300_PFMET100_);
      tree_->SetBranchStatus("HLT_PFHT300_PFMET110", 1); tree_->SetBranchAddress("HLT_PFHT300_PFMET110", &HLT_PFHT300_PFMET110_);
      tree_->SetBranchStatus("HLT_PFHT550_4JetPt50", 1); tree_->SetBranchAddress("HLT_PFHT550_4JetPt50", &HLT_PFHT550_4JetPt50_);
      tree_->SetBranchStatus("HLT_PFHT650_4JetPt50", 1); tree_->SetBranchAddress("HLT_PFHT650_4JetPt50", &HLT_PFHT650_4JetPt50_);
      tree_->SetBranchStatus("HLT_PFHT750_4JetPt50", 1); tree_->SetBranchAddress("HLT_PFHT750_4JetPt50", &HLT_PFHT750_4JetPt50_);
      tree_->SetBranchStatus("HLT_PFHT750_4JetPt70", 1); tree_->SetBranchAddress("HLT_PFHT750_4JetPt70", &HLT_PFHT750_4JetPt70_);
      tree_->SetBranchStatus("HLT_PFHT750_4JetPt80", 1); tree_->SetBranchAddress("HLT_PFHT750_4JetPt80", &HLT_PFHT750_4JetPt80_);
      tree_->SetBranchStatus("HLT_PFHT800_4JetPt50", 1); tree_->SetBranchAddress("HLT_PFHT800_4JetPt50", &HLT_PFHT800_4JetPt50_);
      tree_->SetBranchStatus("HLT_PFHT850_4JetPt50", 1); tree_->SetBranchAddress("HLT_PFHT850_4JetPt50", &HLT_PFHT850_4JetPt50_);
      tree_->SetBranchStatus("HLT_PFJet15_NoCaloMatched", 1); tree_->SetBranchAddress("HLT_PFJet15_NoCaloMatched", &HLT_PFJet15_NoCaloMatched_);
      tree_->SetBranchStatus("HLT_PFJet25_NoCaloMatched", 1); tree_->SetBranchAddress("HLT_PFJet25_NoCaloMatched", &HLT_PFJet25_NoCaloMatched_);
      tree_->SetBranchStatus("HLT_DiPFJet15_NoCaloMatched", 1); tree_->SetBranchAddress("HLT_DiPFJet15_NoCaloMatched", &HLT_DiPFJet15_NoCaloMatched_);
      tree_->SetBranchStatus("HLT_DiPFJet25_NoCaloMatched", 1); tree_->SetBranchAddress("HLT_DiPFJet25_NoCaloMatched", &HLT_DiPFJet25_NoCaloMatched_);
      tree_->SetBranchStatus("HLT_DiPFJet15_FBEta3_NoCaloMatched", 1); tree_->SetBranchAddress("HLT_DiPFJet15_FBEta3_NoCaloMatched", &HLT_DiPFJet15_FBEta3_NoCaloMatched_);
      tree_->SetBranchStatus("HLT_DiPFJet25_FBEta3_NoCaloMatched", 1); tree_->SetBranchAddress("HLT_DiPFJet25_FBEta3_NoCaloMatched", &HLT_DiPFJet25_FBEta3_NoCaloMatched_);
      tree_->SetBranchStatus("HLT_DiPFJetAve15_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve15_HFJEC", &HLT_DiPFJetAve15_HFJEC_);
      tree_->SetBranchStatus("HLT_DiPFJetAve25_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve25_HFJEC", &HLT_DiPFJetAve25_HFJEC_);
      tree_->SetBranchStatus("HLT_DiPFJetAve35_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve35_HFJEC", &HLT_DiPFJetAve35_HFJEC_);
      tree_->SetBranchStatus("HLT_AK8PFJet40", 1); tree_->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40_);
      tree_->SetBranchStatus("HLT_AK8PFJet60", 1); tree_->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60_);
      tree_->SetBranchStatus("HLT_AK8PFJet80", 1); tree_->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80_);
      tree_->SetBranchStatus("HLT_AK8PFJet140", 1); tree_->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140_);
      tree_->SetBranchStatus("HLT_AK8PFJet200", 1); tree_->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200_);
      tree_->SetBranchStatus("HLT_AK8PFJet260", 1); tree_->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260_);
      tree_->SetBranchStatus("HLT_AK8PFJet320", 1); tree_->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320_);
      tree_->SetBranchStatus("HLT_AK8PFJet400", 1); tree_->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400_);
      tree_->SetBranchStatus("HLT_AK8PFJet450", 1); tree_->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450_);
      tree_->SetBranchStatus("HLT_AK8PFJet500", 1); tree_->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500_);
      tree_->SetBranchStatus("HLT_PFJet40", 1); tree_->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40_);
      tree_->SetBranchStatus("HLT_PFJet60", 1); tree_->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60_);
      tree_->SetBranchStatus("HLT_PFJet80", 1); tree_->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80_);
      tree_->SetBranchStatus("HLT_PFJet140", 1); tree_->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140_);
      tree_->SetBranchStatus("HLT_PFJet200", 1); tree_->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200_);
      tree_->SetBranchStatus("HLT_PFJet260", 1); tree_->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260_);
      tree_->SetBranchStatus("HLT_PFJet320", 1); tree_->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320_);
      tree_->SetBranchStatus("HLT_PFJet400", 1); tree_->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400_);
      tree_->SetBranchStatus("HLT_PFJet450", 1); tree_->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450_);
      tree_->SetBranchStatus("HLT_PFJet500", 1); tree_->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500_);
      tree_->SetBranchStatus("HLT_DiPFJetAve40", 1); tree_->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40_);
      tree_->SetBranchStatus("HLT_DiPFJetAve60", 1); tree_->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60_);
      tree_->SetBranchStatus("HLT_DiPFJetAve80", 1); tree_->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80_);
      tree_->SetBranchStatus("HLT_DiPFJetAve140", 1); tree_->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140_);
      tree_->SetBranchStatus("HLT_DiPFJetAve200", 1); tree_->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200_);
      tree_->SetBranchStatus("HLT_DiPFJetAve260", 1); tree_->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260_);
      tree_->SetBranchStatus("HLT_DiPFJetAve320", 1); tree_->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320_);
      tree_->SetBranchStatus("HLT_DiPFJetAve400", 1); tree_->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400_);
      tree_->SetBranchStatus("HLT_DiPFJetAve500", 1); tree_->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500_);
      tree_->SetBranchStatus("HLT_DiPFJetAve60_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC_);
      tree_->SetBranchStatus("HLT_DiPFJetAve80_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC_);
      tree_->SetBranchStatus("HLT_DiPFJetAve100_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC_);
      tree_->SetBranchStatus("HLT_DiPFJetAve160_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC_);
      tree_->SetBranchStatus("HLT_DiPFJetAve220_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC_);
      tree_->SetBranchStatus("HLT_DiPFJetAve300_HFJEC", 1); tree_->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC_);
      tree_->SetBranchStatus("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140", 1); tree_->SetBranchAddress("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140", &HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140_);
      tree_->SetBranchStatus("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80", 1); tree_->SetBranchAddress("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80", &HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80_);
      tree_->SetBranchStatus("HLT_DiCentralPFJet170", 1); tree_->SetBranchAddress("HLT_DiCentralPFJet170", &HLT_DiCentralPFJet170_);
      tree_->SetBranchStatus("HLT_SingleCentralPFJet170_CFMax0p1", 1); tree_->SetBranchAddress("HLT_SingleCentralPFJet170_CFMax0p1", &HLT_SingleCentralPFJet170_CFMax0p1_);
      tree_->SetBranchStatus("HLT_DiCentralPFJet170_CFMax0p1", 1); tree_->SetBranchAddress("HLT_DiCentralPFJet170_CFMax0p1", &HLT_DiCentralPFJet170_CFMax0p1_);
      tree_->SetBranchStatus("HLT_DiCentralPFJet220_CFMax0p3", 1); tree_->SetBranchAddress("HLT_DiCentralPFJet220_CFMax0p3", &HLT_DiCentralPFJet220_CFMax0p3_);
      tree_->SetBranchStatus("HLT_DiCentralPFJet330_CFMax0p5", 1); tree_->SetBranchAddress("HLT_DiCentralPFJet330_CFMax0p5", &HLT_DiCentralPFJet330_CFMax0p5_);
      tree_->SetBranchStatus("HLT_DiCentralPFJet430", 1); tree_->SetBranchAddress("HLT_DiCentralPFJet430", &HLT_DiCentralPFJet430_);
      tree_->SetBranchStatus("HLT_PFHT125", 1); tree_->SetBranchAddress("HLT_PFHT125", &HLT_PFHT125_);
      tree_->SetBranchStatus("HLT_PFHT200", 1); tree_->SetBranchAddress("HLT_PFHT200", &HLT_PFHT200_);
      tree_->SetBranchStatus("HLT_PFHT250", 1); tree_->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250_);
      tree_->SetBranchStatus("HLT_PFHT300", 1); tree_->SetBranchAddress("HLT_PFHT300", &HLT_PFHT300_);
      tree_->SetBranchStatus("HLT_PFHT350", 1); tree_->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350_);
      tree_->SetBranchStatus("HLT_PFHT400", 1); tree_->SetBranchAddress("HLT_PFHT400", &HLT_PFHT400_);
      tree_->SetBranchStatus("HLT_PFHT475", 1); tree_->SetBranchAddress("HLT_PFHT475", &HLT_PFHT475_);
      tree_->SetBranchStatus("HLT_PFHT600", 1); tree_->SetBranchAddress("HLT_PFHT600", &HLT_PFHT600_);
      tree_->SetBranchStatus("HLT_PFHT650", 1); tree_->SetBranchAddress("HLT_PFHT650", &HLT_PFHT650_);
      tree_->SetBranchStatus("HLT_PFHT800", 1); tree_->SetBranchAddress("HLT_PFHT800", &HLT_PFHT800_);
      tree_->SetBranchStatus("HLT_PFHT900", 1); tree_->SetBranchAddress("HLT_PFHT900", &HLT_PFHT900_);
      tree_->SetBranchStatus("HLT_PFHT200_PFAlphaT0p51", 1); tree_->SetBranchAddress("HLT_PFHT200_PFAlphaT0p51", &HLT_PFHT200_PFAlphaT0p51_);
      tree_->SetBranchStatus("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57", 1); tree_->SetBranchAddress("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57", &HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57_);
      tree_->SetBranchStatus("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63", 1); tree_->SetBranchAddress("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63", &HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_);
      tree_->SetBranchStatus("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55", 1); tree_->SetBranchAddress("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55", &HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55_);
      tree_->SetBranchStatus("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58", 1); tree_->SetBranchAddress("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58", &HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_);
      tree_->SetBranchStatus("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53", 1); tree_->SetBranchAddress("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53", &HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53_);
      tree_->SetBranchStatus("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54", 1); tree_->SetBranchAddress("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54", &HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_);
      tree_->SetBranchStatus("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52", 1); tree_->SetBranchAddress("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52", &HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52_);
      tree_->SetBranchStatus("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53", 1); tree_->SetBranchAddress("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53", &HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53_);
      tree_->SetBranchStatus("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51", 1); tree_->SetBranchAddress("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51", &HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51_);
      tree_->SetBranchStatus("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52", 1); tree_->SetBranchAddress("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52", &HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52_);
      tree_->SetBranchStatus("HLT_MET60_IsoTrk35_Loose", 1); tree_->SetBranchAddress("HLT_MET60_IsoTrk35_Loose", &HLT_MET60_IsoTrk35_Loose_);
      tree_->SetBranchStatus("HLT_MET75_IsoTrk50", 1); tree_->SetBranchAddress("HLT_MET75_IsoTrk50", &HLT_MET75_IsoTrk50_);
      tree_->SetBranchStatus("HLT_MET90_IsoTrk50", 1); tree_->SetBranchAddress("HLT_MET90_IsoTrk50", &HLT_MET90_IsoTrk50_);
      tree_->SetBranchStatus("HLT_PFMET120_BTagCSV_p067", 1); tree_->SetBranchAddress("HLT_PFMET120_BTagCSV_p067", &HLT_PFMET120_BTagCSV_p067_);
      tree_->SetBranchStatus("HLT_PFMET120_Mu5", 1); tree_->SetBranchAddress("HLT_PFMET120_Mu5", &HLT_PFMET120_Mu5_);
      tree_->SetBranchStatus("HLT_PFMET170_NotCleaned", 1); tree_->SetBranchAddress("HLT_PFMET170_NotCleaned", &HLT_PFMET170_NotCleaned_);
      tree_->SetBranchStatus("HLT_PFMET170_NoiseCleaned", 1); tree_->SetBranchAddress("HLT_PFMET170_NoiseCleaned", &HLT_PFMET170_NoiseCleaned_);
      tree_->SetBranchStatus("HLT_PFMET170_HBHECleaned", 1); tree_->SetBranchAddress("HLT_PFMET170_HBHECleaned", &HLT_PFMET170_HBHECleaned_);
      tree_->SetBranchStatus("HLT_PFMET170_JetIdCleaned", 1); tree_->SetBranchAddress("HLT_PFMET170_JetIdCleaned", &HLT_PFMET170_JetIdCleaned_);
      tree_->SetBranchStatus("HLT_PFMET170_BeamHaloCleaned", 1); tree_->SetBranchAddress("HLT_PFMET170_BeamHaloCleaned", &HLT_PFMET170_BeamHaloCleaned_);
      tree_->SetBranchStatus("HLT_PFMET170_HBHE_BeamHaloCleaned", 1); tree_->SetBranchAddress("HLT_PFMET170_HBHE_BeamHaloCleaned", &HLT_PFMET170_HBHE_BeamHaloCleaned_);
      tree_->SetBranchStatus("HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned", 1); tree_->SetBranchAddress("HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned_);
      tree_->SetBranchStatus("HLT_PFMET90_PFMHT90_IDTight", 1); tree_->SetBranchAddress("HLT_PFMET90_PFMHT90_IDTight", &HLT_PFMET90_PFMHT90_IDTight_);
      tree_->SetBranchStatus("HLT_PFMET100_PFMHT100_IDTight", 1); tree_->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight", &HLT_PFMET100_PFMHT100_IDTight_);
      tree_->SetBranchStatus("HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned", 1); tree_->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned", &HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned_);
      tree_->SetBranchStatus("HLT_PFMET110_PFMHT110_IDTight", 1); tree_->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight_);
      tree_->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight", 1); tree_->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight_);
      tree_->SetBranchStatus("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067", 1); tree_->SetBranchAddress("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067", &HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_);
      tree_->SetBranchStatus("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight", 1); tree_->SetBranchAddress("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight", &HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_);
      tree_->SetBranchStatus("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200", 1); tree_->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200", &HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200_);
      tree_->SetBranchStatus("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460", 1); tree_->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460", &HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460_);
      tree_->SetBranchStatus("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240", 1); tree_->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240", &HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_);
      tree_->SetBranchStatus("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500", 1); tree_->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500", &HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_);
      tree_->SetBranchStatus("HLT_QuadPFJet_VBF", 1); tree_->SetBranchAddress("HLT_QuadPFJet_VBF", &HLT_QuadPFJet_VBF_);
      tree_->SetBranchStatus("HLT_L1_TripleJet_VBF", 1); tree_->SetBranchAddress("HLT_L1_TripleJet_VBF", &HLT_L1_TripleJet_VBF_);
      tree_->SetBranchStatus("HLT_QuadJet45_TripleBTagCSV_p087", 1); tree_->SetBranchAddress("HLT_QuadJet45_TripleBTagCSV_p087", &HLT_QuadJet45_TripleBTagCSV_p087_);
      tree_->SetBranchStatus("HLT_QuadJet45_DoubleBTagCSV_p087", 1); tree_->SetBranchAddress("HLT_QuadJet45_DoubleBTagCSV_p087", &HLT_QuadJet45_DoubleBTagCSV_p087_);
      tree_->SetBranchStatus("HLT_DoubleJet90_Double30_TripleBTagCSV_p087", 1); tree_->SetBranchAddress("HLT_DoubleJet90_Double30_TripleBTagCSV_p087", &HLT_DoubleJet90_Double30_TripleBTagCSV_p087_);
      tree_->SetBranchStatus("HLT_DoubleJet90_Double30_DoubleBTagCSV_p087", 1); tree_->SetBranchAddress("HLT_DoubleJet90_Double30_DoubleBTagCSV_p087", &HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_);
      tree_->SetBranchStatus("HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160", 1); tree_->SetBranchAddress("HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160", &HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_);
      tree_->SetBranchStatus("HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6", 1); tree_->SetBranchAddress("HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6", &HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_);
      tree_->SetBranchStatus("HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172", 1); tree_->SetBranchAddress("HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172", &HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172_);
      tree_->SetBranchStatus("HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6", 1); tree_->SetBranchAddress("HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6", &HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_);
      tree_->SetBranchStatus("HLT_DoubleJetsC100_SingleBTagCSV_p026", 1); tree_->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p026", &HLT_DoubleJetsC100_SingleBTagCSV_p026_);
      tree_->SetBranchStatus("HLT_DoubleJetsC100_SingleBTagCSV_p014", 1); tree_->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p014", &HLT_DoubleJetsC100_SingleBTagCSV_p014_);
      tree_->SetBranchStatus("HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350", 1); tree_->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350", &HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350_);
      tree_->SetBranchStatus("HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350", 1); tree_->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350", &HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350_);
      tree_->SetBranchStatus("HLT_Photon135_PFMET100", 1); tree_->SetBranchAddress("HLT_Photon135_PFMET100", &HLT_Photon135_PFMET100_);
      tree_->SetBranchStatus("HLT_Photon20_CaloIdVL_IsoL", 1); tree_->SetBranchAddress("HLT_Photon20_CaloIdVL_IsoL", &HLT_Photon20_CaloIdVL_IsoL_);
      tree_->SetBranchStatus("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40", 1); tree_->SetBranchAddress("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
      tree_->SetBranchStatus("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF", 1); tree_->SetBranchAddress("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_);
      tree_->SetBranchStatus("HLT_Photon250_NoHE", 1); tree_->SetBranchAddress("HLT_Photon250_NoHE", &HLT_Photon250_NoHE_);
      tree_->SetBranchStatus("HLT_Photon300_NoHE", 1); tree_->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE_);
      tree_->SetBranchStatus("HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60", 1); tree_->SetBranchAddress("HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60", &HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60_);
      tree_->SetBranchStatus("HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15", 1); tree_->SetBranchAddress("HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15", &HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15_);
      tree_->SetBranchStatus("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40", 1); tree_->SetBranchAddress("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
      tree_->SetBranchStatus("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF", 1); tree_->SetBranchAddress("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_);
      tree_->SetBranchStatus("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40", 1); tree_->SetBranchAddress("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
      tree_->SetBranchStatus("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF", 1); tree_->SetBranchAddress("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_);
      tree_->SetBranchStatus("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40", 1); tree_->SetBranchAddress("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
      tree_->SetBranchStatus("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF", 1); tree_->SetBranchAddress("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_);
      tree_->SetBranchStatus("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40", 1); tree_->SetBranchAddress("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
      tree_->SetBranchStatus("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF", 1); tree_->SetBranchAddress("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_);
      tree_->SetBranchStatus("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40", 1); tree_->SetBranchAddress("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
      tree_->SetBranchStatus("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF", 1); tree_->SetBranchAddress("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_);
      tree_->SetBranchStatus("HLT_BTagMu_DiJet20_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_DiJet20_Mu5", &HLT_BTagMu_DiJet20_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_DiJet40_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_DiJet40_Mu5", &HLT_BTagMu_DiJet40_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_DiJet70_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_DiJet70_Mu5", &HLT_BTagMu_DiJet70_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_DiJet110_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_DiJet110_Mu5", &HLT_BTagMu_DiJet110_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_DiJet170_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_DiJet170_Mu5", &HLT_BTagMu_DiJet170_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_Jet300_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_Jet300_Mu5", &HLT_BTagMu_Jet300_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK8Jet300_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5_);
      tree_->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded", 1); tree_->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_);
      tree_->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", 1); tree_->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL", &HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL", &HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL", &HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL", &HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", 1); tree_->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_);
      tree_->SetBranchStatus("HLT_Mu12_Photon25_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL", &HLT_Mu12_Photon25_CaloIdL_);
      tree_->SetBranchStatus("HLT_Mu12_Photon25_CaloIdL_L1ISO", 1); tree_->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL_L1ISO", &HLT_Mu12_Photon25_CaloIdL_L1ISO_);
      tree_->SetBranchStatus("HLT_Mu12_Photon25_CaloIdL_L1OR", 1); tree_->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL_L1OR", &HLT_Mu12_Photon25_CaloIdL_L1OR_);
      tree_->SetBranchStatus("HLT_Mu17_Photon22_CaloIdL_L1ISO", 1); tree_->SetBranchAddress("HLT_Mu17_Photon22_CaloIdL_L1ISO", &HLT_Mu17_Photon22_CaloIdL_L1ISO_);
      tree_->SetBranchStatus("HLT_Mu17_Photon30_CaloIdL_L1ISO", 1); tree_->SetBranchAddress("HLT_Mu17_Photon30_CaloIdL_L1ISO", &HLT_Mu17_Photon30_CaloIdL_L1ISO_);
      tree_->SetBranchStatus("HLT_Mu17_Photon35_CaloIdL_L1ISO", 1); tree_->SetBranchAddress("HLT_Mu17_Photon35_CaloIdL_L1ISO", &HLT_Mu17_Photon35_CaloIdL_L1ISO_);
      tree_->SetBranchStatus("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", 1); tree_->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_);
      tree_->SetBranchStatus("HLT_TripleMu_5_3_3", 1); tree_->SetBranchAddress("HLT_TripleMu_5_3_3", &HLT_TripleMu_5_3_3_);
      tree_->SetBranchStatus("HLT_TripleMu_12_10_5", 1); tree_->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5_);
      tree_->SetBranchStatus("HLT_Mu3er_PFHT140_PFMET125", 1); tree_->SetBranchAddress("HLT_Mu3er_PFHT140_PFMET125", &HLT_Mu3er_PFHT140_PFMET125_);
      tree_->SetBranchStatus("HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067", 1); tree_->SetBranchAddress("HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067", &HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067_);
      tree_->SetBranchStatus("HLT_Mu6_PFHT200_PFMET100", 1); tree_->SetBranchAddress("HLT_Mu6_PFHT200_PFMET100", &HLT_Mu6_PFHT200_PFMET100_);
      tree_->SetBranchStatus("HLT_Mu14er_PFMET100", 1); tree_->SetBranchAddress("HLT_Mu14er_PFMET100", &HLT_Mu14er_PFMET100_);
      tree_->SetBranchStatus("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Ele12_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Ele17_CaloIdL_GsfTrkIdVL", 1); tree_->SetBranchAddress("HLT_Ele17_CaloIdL_GsfTrkIdVL", &HLT_Ele17_CaloIdL_GsfTrkIdVL_);
      tree_->SetBranchStatus("HLT_Ele17_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Ele23_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5", 1); tree_->SetBranchAddress("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5", &HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_);
      tree_->SetBranchStatus("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5", 1); tree_->SetBranchAddress("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5", &HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_);
      tree_->SetBranchStatus("HLT_Photon22", 1); tree_->SetBranchAddress("HLT_Photon22", &HLT_Photon22_);
      tree_->SetBranchStatus("HLT_Photon30", 1); tree_->SetBranchAddress("HLT_Photon30", &HLT_Photon30_);
      tree_->SetBranchStatus("HLT_Photon36", 1); tree_->SetBranchAddress("HLT_Photon36", &HLT_Photon36_);
      tree_->SetBranchStatus("HLT_Photon50", 1); tree_->SetBranchAddress("HLT_Photon50", &HLT_Photon50_);
      tree_->SetBranchStatus("HLT_Photon75", 1); tree_->SetBranchAddress("HLT_Photon75", &HLT_Photon75_);
      tree_->SetBranchStatus("HLT_Photon90", 1); tree_->SetBranchAddress("HLT_Photon90", &HLT_Photon90_);
      tree_->SetBranchStatus("HLT_Photon120", 1); tree_->SetBranchAddress("HLT_Photon120", &HLT_Photon120_);
      tree_->SetBranchStatus("HLT_Photon175", 1); tree_->SetBranchAddress("HLT_Photon175", &HLT_Photon175_);
      tree_->SetBranchStatus("HLT_Photon165_HE10", 1); tree_->SetBranchAddress("HLT_Photon165_HE10", &HLT_Photon165_HE10_);
      tree_->SetBranchStatus("HLT_Photon22_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon22_R9Id90_HE10_IsoM", &HLT_Photon22_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon30_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon30_R9Id90_HE10_IsoM", &HLT_Photon30_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon36_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon36_R9Id90_HE10_IsoM", &HLT_Photon36_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon50_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon90_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon120_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon165_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", 1); tree_->SetBranchAddress("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_);
      tree_->SetBranchStatus("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70", 1); tree_->SetBranchAddress("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70", &HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_);
      tree_->SetBranchStatus("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", 1); tree_->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_);
      tree_->SetBranchStatus("HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55", 1); tree_->SetBranchAddress("HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55", &HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55_);
      tree_->SetBranchStatus("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", 1); tree_->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_);
      tree_->SetBranchStatus("HLT_Dimuon0_Jpsi_Muon", 1); tree_->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon", &HLT_Dimuon0_Jpsi_Muon_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_Muon", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon", &HLT_Dimuon0_Upsilon_Muon_);
      tree_->SetBranchStatus("HLT_QuadMuon0_Dimuon0_Jpsi", 1); tree_->SetBranchAddress("HLT_QuadMuon0_Dimuon0_Jpsi", &HLT_QuadMuon0_Dimuon0_Jpsi_);
      tree_->SetBranchStatus("HLT_QuadMuon0_Dimuon0_Upsilon", 1); tree_->SetBranchAddress("HLT_QuadMuon0_Dimuon0_Upsilon", &HLT_QuadMuon0_Dimuon0_Upsilon_);
      tree_->SetBranchStatus("HLT_Rsq0p25_Calo", 1); tree_->SetBranchAddress("HLT_Rsq0p25_Calo", &HLT_Rsq0p25_Calo_);
      tree_->SetBranchStatus("HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo", 1); tree_->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo", &HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo_);
      tree_->SetBranchStatus("HLT_RsqMR240_Rsq0p09_MR200_Calo", 1); tree_->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_Calo", &HLT_RsqMR240_Rsq0p09_MR200_Calo_);
      tree_->SetBranchStatus("HLT_Rsq0p25", 1); tree_->SetBranchAddress("HLT_Rsq0p25", &HLT_Rsq0p25_);
      tree_->SetBranchStatus("HLT_Rsq0p30", 1); tree_->SetBranchAddress("HLT_Rsq0p30", &HLT_Rsq0p30_);
      tree_->SetBranchStatus("HLT_RsqMR240_Rsq0p09_MR200", 1); tree_->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200", &HLT_RsqMR240_Rsq0p09_MR200_);
      tree_->SetBranchStatus("HLT_RsqMR240_Rsq0p09_MR200_4jet", 1); tree_->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet", &HLT_RsqMR240_Rsq0p09_MR200_4jet_);
      tree_->SetBranchStatus("HLT_RsqMR270_Rsq0p09_MR200", 1); tree_->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200", &HLT_RsqMR270_Rsq0p09_MR200_);
      tree_->SetBranchStatus("HLT_RsqMR270_Rsq0p09_MR200_4jet", 1); tree_->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200_4jet", &HLT_RsqMR270_Rsq0p09_MR200_4jet_);
      tree_->SetBranchStatus("HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200", 1); tree_->SetBranchAddress("HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200", &HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200_);
      tree_->SetBranchStatus("HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", 1); tree_->SetBranchAddress("HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_);
      tree_->SetBranchStatus("HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", 1); tree_->SetBranchAddress("HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_);
      tree_->SetBranchStatus("HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", 1); tree_->SetBranchAddress("HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_);
      tree_->SetBranchStatus("HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", 1); tree_->SetBranchAddress("HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_);
      tree_->SetBranchStatus("HLT_HT200_DisplacedDijet40_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT200_DisplacedDijet40_DisplacedTrack", &HLT_HT200_DisplacedDijet40_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT250_DisplacedDijet40_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT250_DisplacedDijet40_DisplacedTrack", &HLT_HT250_DisplacedDijet40_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT350_DisplacedDijet40_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT350_DisplacedDijet40_DisplacedTrack", &HLT_HT350_DisplacedDijet40_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT350_DisplacedDijet80_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT350_DisplacedDijet80_DisplacedTrack", &HLT_HT350_DisplacedDijet80_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack", &HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT350_DisplacedDijet40_Inclusive", 1); tree_->SetBranchAddress("HLT_HT350_DisplacedDijet40_Inclusive", &HLT_HT350_DisplacedDijet40_Inclusive_);
      tree_->SetBranchStatus("HLT_HT400_DisplacedDijet40_Inclusive", 1); tree_->SetBranchAddress("HLT_HT400_DisplacedDijet40_Inclusive", &HLT_HT400_DisplacedDijet40_Inclusive_);
      tree_->SetBranchStatus("HLT_HT500_DisplacedDijet40_Inclusive", 1); tree_->SetBranchAddress("HLT_HT500_DisplacedDijet40_Inclusive", &HLT_HT500_DisplacedDijet40_Inclusive_);
      tree_->SetBranchStatus("HLT_HT550_DisplacedDijet40_Inclusive", 1); tree_->SetBranchAddress("HLT_HT550_DisplacedDijet40_Inclusive", &HLT_HT550_DisplacedDijet40_Inclusive_);
      tree_->SetBranchStatus("HLT_HT550_DisplacedDijet80_Inclusive", 1); tree_->SetBranchAddress("HLT_HT550_DisplacedDijet80_Inclusive", &HLT_HT550_DisplacedDijet80_Inclusive_);
      tree_->SetBranchStatus("HLT_HT650_DisplacedDijet80_Inclusive", 1); tree_->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive", &HLT_HT650_DisplacedDijet80_Inclusive_);
      tree_->SetBranchStatus("HLT_HT750_DisplacedDijet80_Inclusive", 1); tree_->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive", &HLT_HT750_DisplacedDijet80_Inclusive_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_DisplacedTrack", &HLT_VBF_DisplacedJet40_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5", &HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_TightID_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_TightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_TightID_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_Hadronic", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_Hadronic", &HLT_VBF_DisplacedJet40_Hadronic_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack", &HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_TightID_Hadronic", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_TightID_Hadronic", &HLT_VBF_DisplacedJet40_TightID_Hadronic_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_VTightID_Hadronic", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_VTightID_Hadronic", &HLT_VBF_DisplacedJet40_VTightID_Hadronic_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_VVTightID_Hadronic", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_VVTightID_Hadronic", &HLT_VBF_DisplacedJet40_VVTightID_Hadronic_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_);
      tree_->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight", 1); tree_->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight_);
      tree_->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight", 1); tree_->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight_);
      tree_->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", 1); tree_->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_);
      tree_->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", 1); tree_->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_);
      tree_->SetBranchStatus("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200", 1); tree_->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200", &HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_);
      tree_->SetBranchStatus("HLT_Photon90_CaloIdL_PFHT500", 1); tree_->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT500", &HLT_Photon90_CaloIdL_PFHT500_);
      tree_->SetBranchStatus("HLT_DoubleMu8_Mass8_PFHT250", 1); tree_->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT250", &HLT_DoubleMu8_Mass8_PFHT250_);
      tree_->SetBranchStatus("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250", 1); tree_->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250_);
      tree_->SetBranchStatus("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250", 1); tree_->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250_);
      tree_->SetBranchStatus("HLT_DoubleMu8_Mass8_PFHT300", 1); tree_->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT300", &HLT_DoubleMu8_Mass8_PFHT300_);
      tree_->SetBranchStatus("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300", 1); tree_->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_);
      tree_->SetBranchStatus("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300", 1); tree_->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_);
      tree_->SetBranchStatus("HLT_Mu10_CentralPFJet30_BTagCSV_p13", 1); tree_->SetBranchAddress("HLT_Mu10_CentralPFJet30_BTagCSV_p13", &HLT_Mu10_CentralPFJet30_BTagCSV_p13_);
      tree_->SetBranchStatus("HLT_DoubleMu3_PFMET50", 1); tree_->SetBranchAddress("HLT_DoubleMu3_PFMET50", &HLT_DoubleMu3_PFMET50_);
      tree_->SetBranchStatus("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13", 1); tree_->SetBranchAddress("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13", &HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400", &HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT350_PFMET50", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350_PFMET50", &HLT_Ele15_IsoVVVL_PFHT350_PFMET50_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT600", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT350", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350", &HLT_Ele15_IsoVVVL_PFHT350_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT400_PFMET50", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT400_PFMET50", &HLT_Ele15_IsoVVVL_PFHT400_PFMET50_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT400", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT400", &HLT_Ele15_IsoVVVL_PFHT400_);
      tree_->SetBranchStatus("HLT_Ele50_IsoVVVL_PFHT400", 1); tree_->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT400", &HLT_Ele50_IsoVVVL_PFHT400_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_);
      tree_->SetBranchStatus("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", 1); tree_->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400", &HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT350_PFMET50", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350_PFMET50", &HLT_Mu15_IsoVVVL_PFHT350_PFMET50_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT600", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT350", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350", &HLT_Mu15_IsoVVVL_PFHT350_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT400_PFMET50", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400_PFMET50", &HLT_Mu15_IsoVVVL_PFHT400_PFMET50_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT400", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400", &HLT_Mu15_IsoVVVL_PFHT400_);
      tree_->SetBranchStatus("HLT_Mu50_IsoVVVL_PFHT400", 1); tree_->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT400", &HLT_Mu50_IsoVVVL_PFHT400_);
      tree_->SetBranchStatus("HLT_Dimuon16_Jpsi", 1); tree_->SetBranchAddress("HLT_Dimuon16_Jpsi", &HLT_Dimuon16_Jpsi_);
      tree_->SetBranchStatus("HLT_Dimuon10_Jpsi_Barrel", 1); tree_->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel", &HLT_Dimuon10_Jpsi_Barrel_);
      tree_->SetBranchStatus("HLT_Dimuon8_PsiPrime_Barrel", 1); tree_->SetBranchAddress("HLT_Dimuon8_PsiPrime_Barrel", &HLT_Dimuon8_PsiPrime_Barrel_);
      tree_->SetBranchStatus("HLT_Dimuon8_Upsilon_Barrel", 1); tree_->SetBranchAddress("HLT_Dimuon8_Upsilon_Barrel", &HLT_Dimuon8_Upsilon_Barrel_);
      tree_->SetBranchStatus("HLT_Dimuon0_Phi_Barrel", 1); tree_->SetBranchAddress("HLT_Dimuon0_Phi_Barrel", &HLT_Dimuon0_Phi_Barrel_);
      tree_->SetBranchStatus("HLT_Mu16_TkMu0_dEta18_Onia", 1); tree_->SetBranchAddress("HLT_Mu16_TkMu0_dEta18_Onia", &HLT_Mu16_TkMu0_dEta18_Onia_);
      tree_->SetBranchStatus("HLT_Mu16_TkMu0_dEta18_Phi", 1); tree_->SetBranchAddress("HLT_Mu16_TkMu0_dEta18_Phi", &HLT_Mu16_TkMu0_dEta18_Phi_);
      tree_->SetBranchStatus("HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_Mu8", 1); tree_->SetBranchAddress("HLT_Mu8", &HLT_Mu8_);
      tree_->SetBranchStatus("HLT_Mu17", 1); tree_->SetBranchAddress("HLT_Mu17", &HLT_Mu17_);
      tree_->SetBranchStatus("HLT_Mu3_PFJet40", 1); tree_->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40_);
      tree_->SetBranchStatus("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele12_CaloIdM_TrackIdM_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele12_CaloIdM_TrackIdM_PFJet30", &HLT_Ele12_CaloIdM_TrackIdM_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140", 1); tree_->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_);
      tree_->SetBranchStatus("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", 1); tree_->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_);
      tree_->SetBranchStatus("HLT_PFHT400_SixJet30_DoubleBTagCSV_p056", 1); tree_->SetBranchAddress("HLT_PFHT400_SixJet30_DoubleBTagCSV_p056", &HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_);
      tree_->SetBranchStatus("HLT_PFHT450_SixJet40_BTagCSV_p056", 1); tree_->SetBranchAddress("HLT_PFHT450_SixJet40_BTagCSV_p056", &HLT_PFHT450_SixJet40_BTagCSV_p056_);
      tree_->SetBranchStatus("HLT_PFHT400_SixJet30", 1); tree_->SetBranchAddress("HLT_PFHT400_SixJet30", &HLT_PFHT400_SixJet30_);
      tree_->SetBranchStatus("HLT_PFHT450_SixJet40", 1); tree_->SetBranchAddress("HLT_PFHT450_SixJet40", &HLT_PFHT450_SixJet40_);
      tree_->SetBranchStatus("HLT_Ele115_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_Mu55", 1); tree_->SetBranchAddress("HLT_Mu55", &HLT_Mu55_);
      tree_->SetBranchStatus("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15", 1); tree_->SetBranchAddress("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15", &HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_);
      tree_->SetBranchStatus("HLT_Photon90_CaloIdL_PFHT600", 1); tree_->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT600", &HLT_Photon90_CaloIdL_PFHT600_);
      tree_->SetBranchStatus("HLT_PixelTracks_Multiplicity60ForEndOfFill", 1); tree_->SetBranchAddress("HLT_PixelTracks_Multiplicity60ForEndOfFill", &HLT_PixelTracks_Multiplicity60ForEndOfFill_);
      tree_->SetBranchStatus("HLT_PixelTracks_Multiplicity85ForEndOfFill", 1); tree_->SetBranchAddress("HLT_PixelTracks_Multiplicity85ForEndOfFill", &HLT_PixelTracks_Multiplicity85ForEndOfFill_);
      tree_->SetBranchStatus("HLT_PixelTracks_Multiplicity110ForEndOfFill", 1); tree_->SetBranchAddress("HLT_PixelTracks_Multiplicity110ForEndOfFill", &HLT_PixelTracks_Multiplicity110ForEndOfFill_);
      tree_->SetBranchStatus("HLT_PixelTracks_Multiplicity135ForEndOfFill", 1); tree_->SetBranchAddress("HLT_PixelTracks_Multiplicity135ForEndOfFill", &HLT_PixelTracks_Multiplicity135ForEndOfFill_);
      tree_->SetBranchStatus("HLT_PixelTracks_Multiplicity160ForEndOfFill", 1); tree_->SetBranchAddress("HLT_PixelTracks_Multiplicity160ForEndOfFill", &HLT_PixelTracks_Multiplicity160ForEndOfFill_);
      tree_->SetBranchStatus("HLT_FullTracks_Multiplicity80", 1); tree_->SetBranchAddress("HLT_FullTracks_Multiplicity80", &HLT_FullTracks_Multiplicity80_);
      tree_->SetBranchStatus("HLT_FullTracks_Multiplicity100", 1); tree_->SetBranchAddress("HLT_FullTracks_Multiplicity100", &HLT_FullTracks_Multiplicity100_);
      tree_->SetBranchStatus("HLT_FullTracks_Multiplicity130", 1); tree_->SetBranchAddress("HLT_FullTracks_Multiplicity130", &HLT_FullTracks_Multiplicity130_);
      tree_->SetBranchStatus("HLT_FullTracks_Multiplicity150", 1); tree_->SetBranchAddress("HLT_FullTracks_Multiplicity150", &HLT_FullTracks_Multiplicity150_);
      tree_->SetBranchStatus("HLT_ECALHT800", 1); tree_->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800_);
      tree_->SetBranchStatus("HLT_DiSC30_18_EIso_AND_HE_Mass70", 1); tree_->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70_);
      tree_->SetBranchStatus("HLT_Photon125", 1); tree_->SetBranchAddress("HLT_Photon125", &HLT_Photon125_);
      tree_->SetBranchStatus("HLT_MET100", 1); tree_->SetBranchAddress("HLT_MET100", &HLT_MET100_);
      tree_->SetBranchStatus("HLT_MET150", 1); tree_->SetBranchAddress("HLT_MET150", &HLT_MET150_);
      tree_->SetBranchStatus("HLT_MET200", 1); tree_->SetBranchAddress("HLT_MET200", &HLT_MET200_);
      tree_->SetBranchStatus("HLT_Ele27_HighEta_Ele20_Mass55", 1); tree_->SetBranchAddress("HLT_Ele27_HighEta_Ele20_Mass55", &HLT_Ele27_HighEta_Ele20_Mass55_);
      tree_->SetBranchStatus("HLT_L1FatEvents", 1); tree_->SetBranchAddress("HLT_L1FatEvents", &HLT_L1FatEvents_);
      tree_->SetBranchStatus("HLT_Physics", 1); tree_->SetBranchAddress("HLT_Physics", &HLT_Physics_);
      tree_->SetBranchStatus("HLT_L1FatEvents_part0", 1); tree_->SetBranchAddress("HLT_L1FatEvents_part0", &HLT_L1FatEvents_part0_);
      tree_->SetBranchStatus("HLT_L1FatEvents_part1", 1); tree_->SetBranchAddress("HLT_L1FatEvents_part1", &HLT_L1FatEvents_part1_);
      tree_->SetBranchStatus("HLT_L1FatEvents_part2", 1); tree_->SetBranchAddress("HLT_L1FatEvents_part2", &HLT_L1FatEvents_part2_);
      tree_->SetBranchStatus("HLT_L1FatEvents_part3", 1); tree_->SetBranchAddress("HLT_L1FatEvents_part3", &HLT_L1FatEvents_part3_);
      tree_->SetBranchStatus("HLT_Random", 1); tree_->SetBranchAddress("HLT_Random", &HLT_Random_);
      tree_->SetBranchStatus("HLT_ZeroBias", 1); tree_->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias_);
      tree_->SetBranchStatus("HLT_AK4CaloJet30", 1); tree_->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30_);
      tree_->SetBranchStatus("HLT_AK4CaloJet40", 1); tree_->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40_);
      tree_->SetBranchStatus("HLT_AK4CaloJet50", 1); tree_->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50_);
      tree_->SetBranchStatus("HLT_AK4CaloJet80", 1); tree_->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80_);
      tree_->SetBranchStatus("HLT_AK4CaloJet100", 1); tree_->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100_);
      tree_->SetBranchStatus("HLT_AK4PFJet30", 1); tree_->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30_);
      tree_->SetBranchStatus("HLT_AK4PFJet50", 1); tree_->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50_);
      tree_->SetBranchStatus("HLT_AK4PFJet80", 1); tree_->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80_);
      tree_->SetBranchStatus("HLT_AK4PFJet100", 1); tree_->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100_);
      tree_->SetBranchStatus("HLT_HISinglePhoton10", 1); tree_->SetBranchAddress("HLT_HISinglePhoton10", &HLT_HISinglePhoton10_);
      tree_->SetBranchStatus("HLT_HISinglePhoton15", 1); tree_->SetBranchAddress("HLT_HISinglePhoton15", &HLT_HISinglePhoton15_);
      tree_->SetBranchStatus("HLT_HISinglePhoton20", 1); tree_->SetBranchAddress("HLT_HISinglePhoton20", &HLT_HISinglePhoton20_);
      tree_->SetBranchStatus("HLT_HISinglePhoton40", 1); tree_->SetBranchAddress("HLT_HISinglePhoton40", &HLT_HISinglePhoton40_);
      tree_->SetBranchStatus("HLT_HISinglePhoton60", 1); tree_->SetBranchAddress("HLT_HISinglePhoton60", &HLT_HISinglePhoton60_);
      tree_->SetBranchStatus("HLT_EcalCalibration", 1); tree_->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration_);
      tree_->SetBranchStatus("HLT_HcalCalibration", 1); tree_->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration_);
      tree_->SetBranchStatus("HLT_GlobalRunHPDNoise", 1); tree_->SetBranchAddress("HLT_GlobalRunHPDNoise", &HLT_GlobalRunHPDNoise_);
      tree_->SetBranchStatus("HLT_L1BptxMinus", 1); tree_->SetBranchAddress("HLT_L1BptxMinus", &HLT_L1BptxMinus_);
      tree_->SetBranchStatus("HLT_L1BptxPlus", 1); tree_->SetBranchAddress("HLT_L1BptxPlus", &HLT_L1BptxPlus_);
      tree_->SetBranchStatus("HLT_L1NotBptxOR", 1); tree_->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR_);
      tree_->SetBranchStatus("HLT_L1BeamGasMinus", 1); tree_->SetBranchAddress("HLT_L1BeamGasMinus", &HLT_L1BeamGasMinus_);
      tree_->SetBranchStatus("HLT_L1BeamGasPlus", 1); tree_->SetBranchAddress("HLT_L1BeamGasPlus", &HLT_L1BeamGasPlus_);
      tree_->SetBranchStatus("HLT_L1BptxXOR", 1); tree_->SetBranchAddress("HLT_L1BptxXOR", &HLT_L1BptxXOR_);
      tree_->SetBranchStatus("HLT_L1MinimumBiasHF_OR", 1); tree_->SetBranchAddress("HLT_L1MinimumBiasHF_OR", &HLT_L1MinimumBiasHF_OR_);
      tree_->SetBranchStatus("HLT_L1MinimumBiasHF_AND", 1); tree_->SetBranchAddress("HLT_L1MinimumBiasHF_AND", &HLT_L1MinimumBiasHF_AND_);
      tree_->SetBranchStatus("HLT_HcalNZS", 1); tree_->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS_);
      tree_->SetBranchStatus("HLT_HcalPhiSym", 1); tree_->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym_);
      tree_->SetBranchStatus("HLT_HcalIsolatedbunch", 1); tree_->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch_);
      tree_->SetBranchStatus("HLT_ZeroBias_FirstCollisionAfterAbortGap", 1); tree_->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap_);
      tree_->SetBranchStatus("HLT_ZeroBias_FirstCollisionAfterAbortGap_copy", 1); tree_->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap_copy", &HLT_ZeroBias_FirstCollisionAfterAbortGap_copy_);
      tree_->SetBranchStatus("HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS", 1); tree_->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS", &HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS_);
      tree_->SetBranchStatus("HLT_ZeroBias_IsolatedBunches", 1); tree_->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches_);
      tree_->SetBranchStatus("HLT_ZeroBias_FirstCollisionInTrain", 1); tree_->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain_);
      tree_->SetBranchStatus("HLT_ZeroBias_FirstBXAfterTrain", 1); tree_->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain_);
      tree_->SetBranchStatus("HLT_Photon500", 1); tree_->SetBranchAddress("HLT_Photon500", &HLT_Photon500_);
      tree_->SetBranchStatus("HLT_Photon600", 1); tree_->SetBranchAddress("HLT_Photon600", &HLT_Photon600_);
      tree_->SetBranchStatus("HLT_Mu300", 1); tree_->SetBranchAddress("HLT_Mu300", &HLT_Mu300_);
      tree_->SetBranchStatus("HLT_Mu350", 1); tree_->SetBranchAddress("HLT_Mu350", &HLT_Mu350_);
      tree_->SetBranchStatus("HLT_MET250", 1); tree_->SetBranchAddress("HLT_MET250", &HLT_MET250_);
      tree_->SetBranchStatus("HLT_MET300", 1); tree_->SetBranchAddress("HLT_MET300", &HLT_MET300_);
      tree_->SetBranchStatus("HLT_MET600", 1); tree_->SetBranchAddress("HLT_MET600", &HLT_MET600_);
      tree_->SetBranchStatus("HLT_MET700", 1); tree_->SetBranchAddress("HLT_MET700", &HLT_MET700_);
      tree_->SetBranchStatus("HLT_PFMET300", 1); tree_->SetBranchAddress("HLT_PFMET300", &HLT_PFMET300_);
      tree_->SetBranchStatus("HLT_PFMET400", 1); tree_->SetBranchAddress("HLT_PFMET400", &HLT_PFMET400_);
      tree_->SetBranchStatus("HLT_PFMET500", 1); tree_->SetBranchAddress("HLT_PFMET500", &HLT_PFMET500_);
      tree_->SetBranchStatus("HLT_PFMET600", 1); tree_->SetBranchAddress("HLT_PFMET600", &HLT_PFMET600_);
      tree_->SetBranchStatus("HLT_Ele250_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele250_CaloIdVT_GsfTrkIdT", &HLT_Ele250_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_Ele300_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele300_CaloIdVT_GsfTrkIdT", &HLT_Ele300_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_HT2000", 1); tree_->SetBranchAddress("HLT_HT2000", &HLT_HT2000_);
      tree_->SetBranchStatus("HLT_HT2500", 1); tree_->SetBranchAddress("HLT_HT2500", &HLT_HT2500_);
      tree_->SetBranchStatus("HLT_IsoTrackHE", 1); tree_->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE_);
      tree_->SetBranchStatus("HLT_IsoTrackHB", 1); tree_->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB_);
      tree_->SetBranchStatus("HLTriggerFinalPath", 1); tree_->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath);
      are_HLT_loaded_ = true;
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
  
  void loadMet(){
    if(!are_MET_loaded_){
      tree_->SetBranchStatus("MET_MetUnclustEnUpDeltaX", 1); tree_->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX_);
      tree_->SetBranchStatus("MET_MetUnclustEnUpDeltaY", 1); tree_->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY_);
      tree_->SetBranchStatus("MET_covXX", 1); tree_->SetBranchAddress("MET_covXX", &MET_covXX_);
      tree_->SetBranchStatus("MET_covXY", 1); tree_->SetBranchAddress("MET_covXY", &MET_covXY_);
      tree_->SetBranchStatus("MET_covYY", 1); tree_->SetBranchAddress("MET_covYY", &MET_covYY_);
      tree_->SetBranchStatus("MET_phi", 1); tree_->SetBranchAddress("MET_phi", &MET_phi_);
      tree_->SetBranchStatus("MET_pt", 1); tree_->SetBranchAddress("MET_pt", &MET_pt_);
      tree_->SetBranchStatus("MET_significance", 1); tree_->SetBranchAddress("MET_significance", &MET_significance_);
      tree_->SetBranchStatus("MET_sumEt", 1); tree_->SetBranchAddress("MET_sumEt", &MET_sumEt_);
      tree_->SetBranchStatus("MET_fiducialGenPhi", 1); tree_->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi_);
      tree_->SetBranchStatus("MET_fiducialGenPt", 1); tree_->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt_);
      are_MET_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenmet(){
    if(!are_GenMET_loaded_){
      tree_->SetBranchStatus("GenMET_phi", 1); tree_->SetBranchAddress("GenMET_phi", &GenMET_phi_);
      tree_->SetBranchStatus("GenMET_pt", 1); tree_->SetBranchAddress("GenMET_pt", &GenMET_pt_);
      are_GenMET_loaded_ = true;
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
  
  void loadFlag(){
    if(!are_Flag_loaded_){
      tree_->SetBranchStatus("Flag_HBHENoiseFilter", 1); tree_->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter_);
      tree_->SetBranchStatus("Flag_HBHENoiseIsoFilter", 1); tree_->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter_);
      tree_->SetBranchStatus("Flag_CSCTightHaloFilter", 1); tree_->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter_);
      tree_->SetBranchStatus("Flag_CSCTightHaloTrkMuUnvetoFilter", 1); tree_->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter_);
      tree_->SetBranchStatus("Flag_CSCTightHalo2015Filter", 1); tree_->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter_);
      tree_->SetBranchStatus("Flag_globalTightHalo2016Filter", 1); tree_->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter_);
      tree_->SetBranchStatus("Flag_globalSuperTightHalo2016Filter", 1); tree_->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter_);
      tree_->SetBranchStatus("Flag_HcalStripHaloFilter", 1); tree_->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter_);
      tree_->SetBranchStatus("Flag_hcalLaserEventFilter", 1); tree_->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter_);
      tree_->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1); tree_->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter_);
      tree_->SetBranchStatus("Flag_EcalDeadCellBoundaryEnergyFilter", 1); tree_->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter_);
      tree_->SetBranchStatus("Flag_goodVertices", 1); tree_->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices_);
      tree_->SetBranchStatus("Flag_eeBadScFilter", 1); tree_->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter_);
      tree_->SetBranchStatus("Flag_ecalLaserCorrFilter", 1); tree_->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter_);
      tree_->SetBranchStatus("Flag_trkPOGFilters", 1); tree_->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters_);
      tree_->SetBranchStatus("Flag_chargedHadronTrackResolutionFilter", 1); tree_->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter_);
      tree_->SetBranchStatus("Flag_muonBadTrackFilter", 1); tree_->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter_);
      tree_->SetBranchStatus("Flag_trkPOG_manystripclus53X", 1); tree_->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X_);
      tree_->SetBranchStatus("Flag_trkPOG_toomanystripclus53X", 1); tree_->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X_);
      tree_->SetBranchStatus("Flag_trkPOG_logErrorTooManyClusters", 1); tree_->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters_);
      tree_->SetBranchStatus("Flag_METFilters", 1); tree_->SetBranchAddress("Flag_METFilters", &Flag_METFilters_);
      are_Flag_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadPileup(){
    if(!are_Pileup_loaded_){
      tree_->SetBranchStatus("Pileup_nTrueInt", 1); tree_->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt_);
      tree_->SetBranchStatus("Pileup_nPU", 1); tree_->SetBranchAddress("Pileup_nPU", &Pileup_nPU_);
      tree_->SetBranchStatus("Pileup_sumEOOT", 1); tree_->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT_);
      tree_->SetBranchStatus("Pileup_sumLOOT", 1); tree_->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT_);
      are_Pileup_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadLheweight(){
    if(!are_LHEWeight_loaded_){
      tree_->SetBranchStatus("LHEWeight_originalXWGTUP", 1); tree_->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP_);
      are_LHEWeight_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  

  const vector<Je>& Jet(){
    if(Jet_.size() > 0) return Jet_;
    loadJes();
  	Jet_.reserve(Jet_area_->size());
    auto it_Jet_area_ = Jet_area_->cbegin();
    auto it_Jet_btagCMVA_ = Jet_btagCMVA_->cbegin();
    auto it_Jet_btagCSVV2_ = Jet_btagCSVV2_->cbegin();
    auto it_Jet_btagDeepB_ = Jet_btagDeepB_->cbegin();
    auto it_Jet_btagDeepC_ = Jet_btagDeepC_->cbegin();
    auto it_Jet_chEmEF_ = Jet_chEmEF_->cbegin();
    auto it_Jet_chHEF_ = Jet_chHEF_->cbegin();
    auto it_Jet_eta_ = Jet_eta_->cbegin();
    auto it_Jet_mass_ = Jet_mass_->cbegin();
    auto it_Jet_neEmEF_ = Jet_neEmEF_->cbegin();
    auto it_Jet_neHEF_ = Jet_neHEF_->cbegin();
    auto it_Jet_phi_ = Jet_phi_->cbegin();
    auto it_Jet_pt_ = Jet_pt_->cbegin();
    auto it_Jet_qgl_ = Jet_qgl_->cbegin();
    auto it_Jet_rawFactor_ = Jet_rawFactor_->cbegin();
    auto it_Jet_bReg_ = Jet_bReg_->cbegin();
    auto it_Jet_electronIdx1_ = Jet_electronIdx1_->cbegin();
    auto it_Jet_electronIdx2_ = Jet_electronIdx2_->cbegin();
    auto it_Jet_jetId_ = Jet_jetId_->cbegin();
    auto it_Jet_muonIdx1_ = Jet_muonIdx1_->cbegin();
    auto it_Jet_muonIdx2_ = Jet_muonIdx2_->cbegin();
    auto it_Jet_nConstituents_ = Jet_nConstituents_->cbegin();
    auto it_Jet_nElectrons_ = Jet_nElectrons_->cbegin();
    auto it_Jet_nMuons_ = Jet_nMuons_->cbegin();
    auto it_Jet_puId_ = Jet_puId_->cbegin();
    auto it_Jet_genJetIdx_ = Jet_genJetIdx_->cbegin();
    auto it_Jet_hadronFlavour_ = Jet_hadronFlavour_->cbegin();
    auto it_Jet_partonFlavour_ = Jet_partonFlavour_->cbegin();
    auto it_Jet_cleanmask_ = Jet_cleanmask_->cbegin();
    for(; it_Jet_area_ != Jet_area_->cend(); ){
      Je obj;
      obj.setarea(*it_Jet_area_);
      obj.setbtagCMVA(*it_Jet_btagCMVA_);
      obj.setbtagCSVV2(*it_Jet_btagCSVV2_);
      obj.setbtagDeepB(*it_Jet_btagDeepB_);
      obj.setbtagDeepC(*it_Jet_btagDeepC_);
      obj.setchEmEF(*it_Jet_chEmEF_);
      obj.setchHEF(*it_Jet_chHEF_);
      obj.setneEmEF(*it_Jet_neEmEF_);
      obj.setneHEF(*it_Jet_neHEF_);
      obj.setqgl(*it_Jet_qgl_);
      obj.setrawFactor(*it_Jet_rawFactor_);
      obj.setbReg(*it_Jet_bReg_);
      obj.setelectronIdx1(*it_Jet_electronIdx1_);
      obj.setelectronIdx2(*it_Jet_electronIdx2_);
      obj.setjetId(*it_Jet_jetId_);
      obj.setmuonIdx1(*it_Jet_muonIdx1_);
      obj.setmuonIdx2(*it_Jet_muonIdx2_);
      obj.setnConstituents(*it_Jet_nConstituents_);
      obj.setnElectrons(*it_Jet_nElectrons_);
      obj.setnMuons(*it_Jet_nMuons_);
      obj.setpuId(*it_Jet_puId_);
      obj.setgenJetIdx(*it_Jet_genJetIdx_);
      obj.sethadronFlavour(*it_Jet_hadronFlavour_);
      obj.setpartonFlavour(*it_Jet_partonFlavour_);
      obj.setcleanmask(*it_Jet_cleanmask_);
      obj.setLotentzVector(*it_Jet_pt_, *it_Jet_eta_, *it_Jet_phi_, *it_Jet_mass_);
      Jet_.push_back( obj );
      ++it_Jet_area_;
      ++it_Jet_btagCMVA_;
      ++it_Jet_btagCSVV2_;
      ++it_Jet_btagDeepB_;
      ++it_Jet_btagDeepC_;
      ++it_Jet_chEmEF_;
      ++it_Jet_chHEF_;
      ++it_Jet_eta_;
      ++it_Jet_mass_;
      ++it_Jet_neEmEF_;
      ++it_Jet_neHEF_;
      ++it_Jet_phi_;
      ++it_Jet_pt_;
      ++it_Jet_qgl_;
      ++it_Jet_rawFactor_;
      ++it_Jet_bReg_;
      ++it_Jet_electronIdx1_;
      ++it_Jet_electronIdx2_;
      ++it_Jet_jetId_;
      ++it_Jet_muonIdx1_;
      ++it_Jet_muonIdx2_;
      ++it_Jet_nConstituents_;
      ++it_Jet_nElectrons_;
      ++it_Jet_nMuons_;
      ++it_Jet_puId_;
      ++it_Jet_genJetIdx_;
      ++it_Jet_hadronFlavour_;
      ++it_Jet_partonFlavour_;
      ++it_Jet_cleanmask_;
    }
    return Jet_;
  }
  
  const vector<Genjetak>& GenJetAK8(){
    if(GenJetAK8_.size() > 0) return GenJetAK8_;
    loadGenjetaks();
  	GenJetAK8_.reserve(GenJetAK8_eta_->size());
    auto it_GenJetAK8_eta_ = GenJetAK8_eta_->cbegin();
    auto it_GenJetAK8_mass_ = GenJetAK8_mass_->cbegin();
    auto it_GenJetAK8_phi_ = GenJetAK8_phi_->cbegin();
    auto it_GenJetAK8_pt_ = GenJetAK8_pt_->cbegin();
    auto it_GenJetAK8_partonFlavour_ = GenJetAK8_partonFlavour_->cbegin();
    auto it_GenJetAK8_hadronFlavour_ = GenJetAK8_hadronFlavour_->cbegin();
    for(; it_GenJetAK8_eta_ != GenJetAK8_eta_->cend(); ){
      Genjetak obj;
      obj.setpartonFlavour(*it_GenJetAK8_partonFlavour_);
      obj.sethadronFlavour(*it_GenJetAK8_hadronFlavour_);
      obj.setLotentzVector(*it_GenJetAK8_pt_, *it_GenJetAK8_eta_, *it_GenJetAK8_phi_, *it_GenJetAK8_mass_);
      GenJetAK8_.push_back( obj );
      ++it_GenJetAK8_eta_;
      ++it_GenJetAK8_mass_;
      ++it_GenJetAK8_phi_;
      ++it_GenJetAK8_pt_;
      ++it_GenJetAK8_partonFlavour_;
      ++it_GenJetAK8_hadronFlavour_;
    }
    return GenJetAK8_;
  }
  
  const vector<Genvista>& GenVisTau(){
    if(GenVisTau_.size() > 0) return GenVisTau_;
    loadGenvistas();
  	GenVisTau_.reserve(GenVisTau_eta_->size());
    auto it_GenVisTau_eta_ = GenVisTau_eta_->cbegin();
    auto it_GenVisTau_mass_ = GenVisTau_mass_->cbegin();
    auto it_GenVisTau_phi_ = GenVisTau_phi_->cbegin();
    auto it_GenVisTau_pt_ = GenVisTau_pt_->cbegin();
    auto it_GenVisTau_charge_ = GenVisTau_charge_->cbegin();
    auto it_GenVisTau_genPartIdxMother_ = GenVisTau_genPartIdxMother_->cbegin();
    auto it_GenVisTau_status_ = GenVisTau_status_->cbegin();
    for(; it_GenVisTau_eta_ != GenVisTau_eta_->cend(); ){
      Genvista obj;
      obj.setcharge(*it_GenVisTau_charge_);
      obj.setgenPartIdxMother(*it_GenVisTau_genPartIdxMother_);
      obj.setstatus(*it_GenVisTau_status_);
      obj.setLotentzVector(*it_GenVisTau_pt_, *it_GenVisTau_eta_, *it_GenVisTau_phi_, *it_GenVisTau_mass_);
      GenVisTau_.push_back( obj );
      ++it_GenVisTau_eta_;
      ++it_GenVisTau_mass_;
      ++it_GenVisTau_phi_;
      ++it_GenVisTau_pt_;
      ++it_GenVisTau_charge_;
      ++it_GenVisTau_genPartIdxMother_;
      ++it_GenVisTau_status_;
    }
    return GenVisTau_;
  }
  
  const Calomet CaloMET(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadCalomet();
  
    Calomet obj;
    obj.setsumEt(CaloMET_sumEt_);
  
    return obj;
  }
  
  const vector<Gendressedlepto>& GenDressedLepton(){
    if(GenDressedLepton_.size() > 0) return GenDressedLepton_;
    loadGendressedleptos();
  	GenDressedLepton_.reserve(GenDressedLepton_eta_->size());
    auto it_GenDressedLepton_eta_ = GenDressedLepton_eta_->cbegin();
    auto it_GenDressedLepton_mass_ = GenDressedLepton_mass_->cbegin();
    auto it_GenDressedLepton_phi_ = GenDressedLepton_phi_->cbegin();
    auto it_GenDressedLepton_pt_ = GenDressedLepton_pt_->cbegin();
    auto it_GenDressedLepton_pdgId_ = GenDressedLepton_pdgId_->cbegin();
    for(; it_GenDressedLepton_eta_ != GenDressedLepton_eta_->cend(); ){
      Gendressedlepto obj;
      obj.setpdgId(*it_GenDressedLepton_pdgId_);
      obj.setLotentzVector(*it_GenDressedLepton_pt_, *it_GenDressedLepton_eta_, *it_GenDressedLepton_phi_, *it_GenDressedLepton_mass_);
      GenDressedLepton_.push_back( obj );
      ++it_GenDressedLepton_eta_;
      ++it_GenDressedLepton_mass_;
      ++it_GenDressedLepton_phi_;
      ++it_GenDressedLepton_pt_;
      ++it_GenDressedLepton_pdgId_;
    }
    return GenDressedLepton_;
  }
  
  const Pv PV(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadPv();
  
    Pv obj;
    obj.setndof(PV_ndof_);
    obj.setx(PV_x_);
    obj.sety(PV_y_);
    obj.setz(PV_z_);
    obj.setchi2(PV_chi2_);
    obj.setscore(PV_score_);
    obj.setnpvs(PV_npvs_);
    obj.setnpvsGood(PV_npvsGood_);
  
    return obj;
  }
  
  const Generator Generator(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGenerator();
  
    Generator obj;
    obj.setbinvar(Generator_binvar_);
    obj.setscalePDF(Generator_scalePDF_);
    obj.setweight(Generator_weight_);
    obj.setx1(Generator_x1_);
    obj.setx2(Generator_x2_);
    obj.setxpdf1(Generator_xpdf1_);
    obj.setxpdf2(Generator_xpdf2_);
    obj.setid1(Generator_id1_);
    obj.setid2(Generator_id2_);
  
    return obj;
  }
  
  const vector<Trigob>& TrigObj(){
    if(TrigObj_.size() > 0) return TrigObj_;
    loadTrigobs();
  	TrigObj_.reserve(TrigObj_pt_->size());
    auto it_TrigObj_pt_ = TrigObj_pt_->cbegin();
    auto it_TrigObj_eta_ = TrigObj_eta_->cbegin();
    auto it_TrigObj_phi_ = TrigObj_phi_->cbegin();
    auto it_TrigObj_l1pt_ = TrigObj_l1pt_->cbegin();
    auto it_TrigObj_l1pt_2_ = TrigObj_l1pt_2_->cbegin();
    auto it_TrigObj_l2pt_ = TrigObj_l2pt_->cbegin();
    auto it_TrigObj_id_ = TrigObj_id_->cbegin();
    auto it_TrigObj_l1iso_ = TrigObj_l1iso_->cbegin();
    auto it_TrigObj_l1charge_ = TrigObj_l1charge_->cbegin();
    auto it_TrigObj_filterBits_ = TrigObj_filterBits_->cbegin();
    for(; it_TrigObj_pt_ != TrigObj_pt_->cend(); ){
      Trigob obj;
      obj.setl1pt(*it_TrigObj_l1pt_);
      obj.setl1pt(*it_TrigObj_l1pt_2_);
      obj.setl2pt(*it_TrigObj_l2pt_);
      obj.setid(*it_TrigObj_id_);
      obj.setl1iso(*it_TrigObj_l1iso_);
      obj.setl1charge(*it_TrigObj_l1charge_);
      obj.setfilterBits(*it_TrigObj_filterBits_);
      obj.setLotentzVector(*it_TrigObj_pt_, *it_TrigObj_eta_, *it_TrigObj_phi_);
      TrigObj_.push_back( obj );
      ++it_TrigObj_pt_;
      ++it_TrigObj_eta_;
      ++it_TrigObj_phi_;
      ++it_TrigObj_l1pt_;
      ++it_TrigObj_l1pt_2_;
      ++it_TrigObj_l2pt_;
      ++it_TrigObj_id_;
      ++it_TrigObj_l1iso_;
      ++it_TrigObj_l1charge_;
      ++it_TrigObj_filterBits_;
    }
    return TrigObj_;
  }
  
  const vector<Photo>& Photon(){
    if(Photon_.size() > 0) return Photon_;
    loadPhotos();
  	Photon_.reserve(Photon_eCorr_->size());
    auto it_Photon_eCorr_ = Photon_eCorr_->cbegin();
    auto it_Photon_energyErr_ = Photon_energyErr_->cbegin();
    auto it_Photon_eta_ = Photon_eta_->cbegin();
    auto it_Photon_hoe_ = Photon_hoe_->cbegin();
    auto it_Photon_mass_ = Photon_mass_->cbegin();
    auto it_Photon_mvaID_ = Photon_mvaID_->cbegin();
    auto it_Photon_pfRelIso03_all_ = Photon_pfRelIso03_all_->cbegin();
    auto it_Photon_pfRelIso03_chg_ = Photon_pfRelIso03_chg_->cbegin();
    auto it_Photon_phi_ = Photon_phi_->cbegin();
    auto it_Photon_pt_ = Photon_pt_->cbegin();
    auto it_Photon_r9_ = Photon_r9_->cbegin();
    auto it_Photon_sieie_ = Photon_sieie_->cbegin();
    auto it_Photon_charge_ = Photon_charge_->cbegin();
    auto it_Photon_cutBased_ = Photon_cutBased_->cbegin();
    auto it_Photon_electronIdx_ = Photon_electronIdx_->cbegin();
    auto it_Photon_jetIdx_ = Photon_jetIdx_->cbegin();
    auto it_Photon_pdgId_ = Photon_pdgId_->cbegin();
    auto it_Photon_vidNestedWPBitmap_ = Photon_vidNestedWPBitmap_->cbegin();
    auto it_Photon_electronVeto_ = Photon_electronVeto_->cbegin();
    auto it_Photon_mvaID_WP80_ = Photon_mvaID_WP80_->cbegin();
    auto it_Photon_mvaID_WP90_ = Photon_mvaID_WP90_->cbegin();
    auto it_Photon_pixelSeed_ = Photon_pixelSeed_->cbegin();
    auto it_Photon_genPartIdx_ = Photon_genPartIdx_->cbegin();
    auto it_Photon_genPartFlav_ = Photon_genPartFlav_->cbegin();
    auto it_Photon_cleanmask_ = Photon_cleanmask_->cbegin();
    for(; it_Photon_eCorr_ != Photon_eCorr_->cend(); ){
      Photo obj;
      obj.seteCorr(*it_Photon_eCorr_);
      obj.setenergyErr(*it_Photon_energyErr_);
      obj.sethoe(*it_Photon_hoe_);
      obj.setmvaID(*it_Photon_mvaID_);
      obj.setpfRelIso03(*it_Photon_pfRelIso03_all_);
      obj.setpfRelIso03(*it_Photon_pfRelIso03_chg_);
      obj.setr9(*it_Photon_r9_);
      obj.setsieie(*it_Photon_sieie_);
      obj.setcharge(*it_Photon_charge_);
      obj.setcutBased(*it_Photon_cutBased_);
      obj.setelectronIdx(*it_Photon_electronIdx_);
      obj.setjetIdx(*it_Photon_jetIdx_);
      obj.setpdgId(*it_Photon_pdgId_);
      obj.setvidNestedWPBitmap(*it_Photon_vidNestedWPBitmap_);
      obj.setelectronVeto(*it_Photon_electronVeto_);
      obj.setmvaID(*it_Photon_mvaID_WP80_);
      obj.setmvaID(*it_Photon_mvaID_WP90_);
      obj.setpixelSeed(*it_Photon_pixelSeed_);
      obj.setgenPartIdx(*it_Photon_genPartIdx_);
      obj.setgenPartFlav(*it_Photon_genPartFlav_);
      obj.setcleanmask(*it_Photon_cleanmask_);
      obj.setLotentzVector(*it_Photon_pt_, *it_Photon_eta_, *it_Photon_phi_, *it_Photon_mass_);
      Photon_.push_back( obj );
      ++it_Photon_eCorr_;
      ++it_Photon_energyErr_;
      ++it_Photon_eta_;
      ++it_Photon_hoe_;
      ++it_Photon_mass_;
      ++it_Photon_mvaID_;
      ++it_Photon_pfRelIso03_all_;
      ++it_Photon_pfRelIso03_chg_;
      ++it_Photon_phi_;
      ++it_Photon_pt_;
      ++it_Photon_r9_;
      ++it_Photon_sieie_;
      ++it_Photon_charge_;
      ++it_Photon_cutBased_;
      ++it_Photon_electronIdx_;
      ++it_Photon_jetIdx_;
      ++it_Photon_pdgId_;
      ++it_Photon_vidNestedWPBitmap_;
      ++it_Photon_electronVeto_;
      ++it_Photon_mvaID_WP80_;
      ++it_Photon_mvaID_WP90_;
      ++it_Photon_pixelSeed_;
      ++it_Photon_genPartIdx_;
      ++it_Photon_genPartFlav_;
      ++it_Photon_cleanmask_;
    }
    return Photon_;
  }
  
  const vector<Genje>& GenJet(){
    if(GenJet_.size() > 0) return GenJet_;
    loadGenjes();
  	GenJet_.reserve(GenJetAK8_eta_->size());
    auto it_GenJetAK8_eta_ = GenJetAK8_eta_->cbegin();
    auto it_GenJetAK8_mass_ = GenJetAK8_mass_->cbegin();
    auto it_GenJetAK8_phi_ = GenJetAK8_phi_->cbegin();
    auto it_GenJetAK8_pt_ = GenJetAK8_pt_->cbegin();
    auto it_GenJet_eta_ = GenJet_eta_->cbegin();
    auto it_GenJet_mass_ = GenJet_mass_->cbegin();
    auto it_GenJet_phi_ = GenJet_phi_->cbegin();
    auto it_GenJet_pt_ = GenJet_pt_->cbegin();
    auto it_GenJetAK8_partonFlavour_ = GenJetAK8_partonFlavour_->cbegin();
    auto it_GenJetAK8_hadronFlavour_ = GenJetAK8_hadronFlavour_->cbegin();
    auto it_GenJet_partonFlavour_ = GenJet_partonFlavour_->cbegin();
    auto it_GenJet_hadronFlavour_ = GenJet_hadronFlavour_->cbegin();
    for(; it_GenJetAK8_eta_ != GenJetAK8_eta_->cend(); ){
      Genje obj;
      obj.setpartonFlavour(*it_GenJetAK8_partonFlavour_);
      obj.sethadronFlavour(*it_GenJetAK8_hadronFlavour_);
      obj.setpartonFlavour(*it_GenJet_partonFlavour_);
      obj.sethadronFlavour(*it_GenJet_hadronFlavour_);
      obj.setLotentzVector(*it_GenJet_pt_, *it_GenJet_eta_, *it_GenJet_phi_, *it_GenJet_mass_);
      GenJet_.push_back( obj );
      ++it_GenJetAK8_eta_;
      ++it_GenJetAK8_mass_;
      ++it_GenJetAK8_phi_;
      ++it_GenJetAK8_pt_;
      ++it_GenJet_eta_;
      ++it_GenJet_mass_;
      ++it_GenJet_phi_;
      ++it_GenJet_pt_;
      ++it_GenJetAK8_partonFlavour_;
      ++it_GenJetAK8_hadronFlavour_;
      ++it_GenJet_partonFlavour_;
      ++it_GenJet_hadronFlavour_;
    }
    return GenJet_;
  }
  
  const Rawmet RawMET(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadRawmet();
  
    Rawmet obj;
    obj.setsumEt(RawMET_sumEt_);
  
    return obj;
  }
  
  const vector<Electro>& Electron(){
    if(Electron_.size() > 0) return Electron_;
    loadElectros();
  	Electron_.reserve(Electron_deltaEtaSC_->size());
    auto it_Electron_deltaEtaSC_ = Electron_deltaEtaSC_->cbegin();
    auto it_Electron_dr03EcalRecHitSumEt_ = Electron_dr03EcalRecHitSumEt_->cbegin();
    auto it_Electron_dr03HcalDepth1TowerSumEt_ = Electron_dr03HcalDepth1TowerSumEt_->cbegin();
    auto it_Electron_dr03TkSumPt_ = Electron_dr03TkSumPt_->cbegin();
    auto it_Electron_dxy_ = Electron_dxy_->cbegin();
    auto it_Electron_dxyErr_ = Electron_dxyErr_->cbegin();
    auto it_Electron_dz_ = Electron_dz_->cbegin();
    auto it_Electron_dzErr_ = Electron_dzErr_->cbegin();
    auto it_Electron_eCorr_ = Electron_eCorr_->cbegin();
    auto it_Electron_eInvMinusPInv_ = Electron_eInvMinusPInv_->cbegin();
    auto it_Electron_energyErr_ = Electron_energyErr_->cbegin();
    auto it_Electron_eta_ = Electron_eta_->cbegin();
    auto it_Electron_hoe_ = Electron_hoe_->cbegin();
    auto it_Electron_ip3d_ = Electron_ip3d_->cbegin();
    auto it_Electron_mass_ = Electron_mass_->cbegin();
    auto it_Electron_miniPFRelIso_all_ = Electron_miniPFRelIso_all_->cbegin();
    auto it_Electron_miniPFRelIso_chg_ = Electron_miniPFRelIso_chg_->cbegin();
    auto it_Electron_mvaSpring16GP_ = Electron_mvaSpring16GP_->cbegin();
    auto it_Electron_mvaSpring16HZZ_ = Electron_mvaSpring16HZZ_->cbegin();
    auto it_Electron_pfRelIso03_all_ = Electron_pfRelIso03_all_->cbegin();
    auto it_Electron_pfRelIso03_chg_ = Electron_pfRelIso03_chg_->cbegin();
    auto it_Electron_phi_ = Electron_phi_->cbegin();
    auto it_Electron_pt_ = Electron_pt_->cbegin();
    auto it_Electron_r9_ = Electron_r9_->cbegin();
    auto it_Electron_sieie_ = Electron_sieie_->cbegin();
    auto it_Electron_sip3d_ = Electron_sip3d_->cbegin();
    auto it_Electron_mvaTTH_ = Electron_mvaTTH_->cbegin();
    auto it_Electron_charge_ = Electron_charge_->cbegin();
    auto it_Electron_cutBased_ = Electron_cutBased_->cbegin();
    auto it_Electron_cutBased_HLTPreSel_ = Electron_cutBased_HLTPreSel_->cbegin();
    auto it_Electron_jetIdx_ = Electron_jetIdx_->cbegin();
    auto it_Electron_pdgId_ = Electron_pdgId_->cbegin();
    auto it_Electron_photonIdx_ = Electron_photonIdx_->cbegin();
    auto it_Electron_tightCharge_ = Electron_tightCharge_->cbegin();
    auto it_Electron_vidNestedWPBitmap_ = Electron_vidNestedWPBitmap_->cbegin();
    auto it_Electron_convVeto_ = Electron_convVeto_->cbegin();
    auto it_Electron_cutBased_HEEP_ = Electron_cutBased_HEEP_->cbegin();
    auto it_Electron_isPFcand_ = Electron_isPFcand_->cbegin();
    auto it_Electron_lostHits_ = Electron_lostHits_->cbegin();
    auto it_Electron_mvaSpring16GP_WP80_ = Electron_mvaSpring16GP_WP80_->cbegin();
    auto it_Electron_mvaSpring16GP_WP90_ = Electron_mvaSpring16GP_WP90_->cbegin();
    auto it_Electron_mvaSpring16HZZ_WPL_ = Electron_mvaSpring16HZZ_WPL_->cbegin();
    auto it_Electron_genPartIdx_ = Electron_genPartIdx_->cbegin();
    auto it_Electron_genPartFlav_ = Electron_genPartFlav_->cbegin();
    auto it_Electron_cleanmask_ = Electron_cleanmask_->cbegin();
    for(; it_Electron_deltaEtaSC_ != Electron_deltaEtaSC_->cend(); ){
      Electro obj;
      obj.setdeltaEtaSC(*it_Electron_deltaEtaSC_);
      obj.setdr03EcalRecHitSumEt(*it_Electron_dr03EcalRecHitSumEt_);
      obj.setdr03HcalDepth1TowerSumEt(*it_Electron_dr03HcalDepth1TowerSumEt_);
      obj.setdr03TkSumPt(*it_Electron_dr03TkSumPt_);
      obj.setdxy(*it_Electron_dxy_);
      obj.setdxyErr(*it_Electron_dxyErr_);
      obj.setdz(*it_Electron_dz_);
      obj.setdzErr(*it_Electron_dzErr_);
      obj.seteCorr(*it_Electron_eCorr_);
      obj.seteInvMinusPInv(*it_Electron_eInvMinusPInv_);
      obj.setenergyErr(*it_Electron_energyErr_);
      obj.sethoe(*it_Electron_hoe_);
      obj.setip3d(*it_Electron_ip3d_);
      obj.setminiPFRelIso(*it_Electron_miniPFRelIso_all_);
      obj.setminiPFRelIso(*it_Electron_miniPFRelIso_chg_);
      obj.setmvaSpring16GP(*it_Electron_mvaSpring16GP_);
      obj.setmvaSpring16HZZ(*it_Electron_mvaSpring16HZZ_);
      obj.setpfRelIso03(*it_Electron_pfRelIso03_all_);
      obj.setpfRelIso03(*it_Electron_pfRelIso03_chg_);
      obj.setr9(*it_Electron_r9_);
      obj.setsieie(*it_Electron_sieie_);
      obj.setsip3d(*it_Electron_sip3d_);
      obj.setmvaTTH(*it_Electron_mvaTTH_);
      obj.setcharge(*it_Electron_charge_);
      obj.setcutBased(*it_Electron_cutBased_);
      obj.setcutBased(*it_Electron_cutBased_HLTPreSel_);
      obj.setjetIdx(*it_Electron_jetIdx_);
      obj.setpdgId(*it_Electron_pdgId_);
      obj.setphotonIdx(*it_Electron_photonIdx_);
      obj.settightCharge(*it_Electron_tightCharge_);
      obj.setvidNestedWPBitmap(*it_Electron_vidNestedWPBitmap_);
      obj.setconvVeto(*it_Electron_convVeto_);
      obj.setcutBased(*it_Electron_cutBased_HEEP_);
      obj.setisPFcand(*it_Electron_isPFcand_);
      obj.setlostHits(*it_Electron_lostHits_);
      obj.setmvaSpring16GP(*it_Electron_mvaSpring16GP_WP80_);
      obj.setmvaSpring16GP(*it_Electron_mvaSpring16GP_WP90_);
      obj.setmvaSpring16HZZ(*it_Electron_mvaSpring16HZZ_WPL_);
      obj.setgenPartIdx(*it_Electron_genPartIdx_);
      obj.setgenPartFlav(*it_Electron_genPartFlav_);
      obj.setcleanmask(*it_Electron_cleanmask_);
      obj.setLotentzVector(*it_Electron_pt_, *it_Electron_eta_, *it_Electron_phi_, *it_Electron_mass_);
      Electron_.push_back( obj );
      ++it_Electron_deltaEtaSC_;
      ++it_Electron_dr03EcalRecHitSumEt_;
      ++it_Electron_dr03HcalDepth1TowerSumEt_;
      ++it_Electron_dr03TkSumPt_;
      ++it_Electron_dxy_;
      ++it_Electron_dxyErr_;
      ++it_Electron_dz_;
      ++it_Electron_dzErr_;
      ++it_Electron_eCorr_;
      ++it_Electron_eInvMinusPInv_;
      ++it_Electron_energyErr_;
      ++it_Electron_eta_;
      ++it_Electron_hoe_;
      ++it_Electron_ip3d_;
      ++it_Electron_mass_;
      ++it_Electron_miniPFRelIso_all_;
      ++it_Electron_miniPFRelIso_chg_;
      ++it_Electron_mvaSpring16GP_;
      ++it_Electron_mvaSpring16HZZ_;
      ++it_Electron_pfRelIso03_all_;
      ++it_Electron_pfRelIso03_chg_;
      ++it_Electron_phi_;
      ++it_Electron_pt_;
      ++it_Electron_r9_;
      ++it_Electron_sieie_;
      ++it_Electron_sip3d_;
      ++it_Electron_mvaTTH_;
      ++it_Electron_charge_;
      ++it_Electron_cutBased_;
      ++it_Electron_cutBased_HLTPreSel_;
      ++it_Electron_jetIdx_;
      ++it_Electron_pdgId_;
      ++it_Electron_photonIdx_;
      ++it_Electron_tightCharge_;
      ++it_Electron_vidNestedWPBitmap_;
      ++it_Electron_convVeto_;
      ++it_Electron_cutBased_HEEP_;
      ++it_Electron_isPFcand_;
      ++it_Electron_lostHits_;
      ++it_Electron_mvaSpring16GP_WP80_;
      ++it_Electron_mvaSpring16GP_WP90_;
      ++it_Electron_mvaSpring16HZZ_WPL_;
      ++it_Electron_genPartIdx_;
      ++it_Electron_genPartFlav_;
      ++it_Electron_cleanmask_;
    }
    return Electron_;
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
  
  const L1Simulation L1simulation(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadL1Simulation();
  
    L1Simulation obj;
    obj.setstep(L1simulation_step_);
  
    return obj;
  }
  
  const vector<Genpar>& GenPart(){
    if(GenPart_.size() > 0) return GenPart_;
    loadGenpars();
  	GenPart_.reserve(GenPart_eta_->size());
    auto it_GenPart_eta_ = GenPart_eta_->cbegin();
    auto it_GenPart_mass_ = GenPart_mass_->cbegin();
    auto it_GenPart_phi_ = GenPart_phi_->cbegin();
    auto it_GenPart_pt_ = GenPart_pt_->cbegin();
    auto it_GenPart_genPartIdxMother_ = GenPart_genPartIdxMother_->cbegin();
    auto it_GenPart_pdgId_ = GenPart_pdgId_->cbegin();
    auto it_GenPart_status_ = GenPart_status_->cbegin();
    auto it_GenPart_statusFlags_ = GenPart_statusFlags_->cbegin();
    for(; it_GenPart_eta_ != GenPart_eta_->cend(); ){
      Genpar obj;
      obj.setgenPartIdxMother(*it_GenPart_genPartIdxMother_);
      obj.setpdgId(*it_GenPart_pdgId_);
      obj.setstatus(*it_GenPart_status_);
      obj.setstatusFlags(*it_GenPart_statusFlags_);
      obj.setLotentzVector(*it_GenPart_pt_, *it_GenPart_eta_, *it_GenPart_phi_, *it_GenPart_mass_);
      GenPart_.push_back( obj );
      ++it_GenPart_eta_;
      ++it_GenPart_mass_;
      ++it_GenPart_phi_;
      ++it_GenPart_pt_;
      ++it_GenPart_genPartIdxMother_;
      ++it_GenPart_pdgId_;
      ++it_GenPart_status_;
      ++it_GenPart_statusFlags_;
    }
    return GenPart_;
  }
  
  const Lhe LHE(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadLhe();
  
    Lhe obj;
    obj.setoriginalXWGTUP(LHEWeight_originalXWGTUP_);
    obj.setLHEPdfWeight(LHEPdfWeight);
    obj.setLHEScaleWeight(LHEScaleWeight);
    obj.setHT(LHE_HT_);
    obj.setHTIncoming(LHE_HTIncoming_);
    obj.setVpt(LHE_Vpt_);
    obj.setNjets(LHE_Njets_);
    obj.setNb(LHE_Nb_);
    obj.setNc(LHE_Nc_);
    obj.setNuds(LHE_Nuds_);
    obj.setNglu(LHE_Nglu_);
    obj.setNpNLO(LHE_NpNLO_);
    obj.setNpLO(LHE_NpLO_);
  
    return obj;
  }
  
  const Tkmet TkMET(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadTkmet();
  
    Tkmet obj;
    obj.setsumEt(TkMET_sumEt_);
  
    return obj;
  }
  
  const vector<Ta>& Tau(){
    if(Tau_.size() > 0) return Tau_;
    loadTas();
  	Tau_.reserve(Tau_chargedIso_->size());
    auto it_Tau_chargedIso_ = Tau_chargedIso_->cbegin();
    auto it_Tau_dxy_ = Tau_dxy_->cbegin();
    auto it_Tau_dz_ = Tau_dz_->cbegin();
    auto it_Tau_eta_ = Tau_eta_->cbegin();
    auto it_Tau_footprintCorr_ = Tau_footprintCorr_->cbegin();
    auto it_Tau_leadTkDeltaEta_ = Tau_leadTkDeltaEta_->cbegin();
    auto it_Tau_leadTkDeltaPhi_ = Tau_leadTkDeltaPhi_->cbegin();
    auto it_Tau_leadTkPtOverTauPt_ = Tau_leadTkPtOverTauPt_->cbegin();
    auto it_Tau_mass_ = Tau_mass_->cbegin();
    auto it_Tau_neutralIso_ = Tau_neutralIso_->cbegin();
    auto it_Tau_phi_ = Tau_phi_->cbegin();
    auto it_Tau_photonsOutsideSignalCone_ = Tau_photonsOutsideSignalCone_->cbegin();
    auto it_Tau_pt_ = Tau_pt_->cbegin();
    auto it_Tau_puCorr_ = Tau_puCorr_->cbegin();
    auto it_Tau_rawAntiEle_ = Tau_rawAntiEle_->cbegin();
    auto it_Tau_rawIso_ = Tau_rawIso_->cbegin();
    auto it_Tau_rawMVAnewDM_ = Tau_rawMVAnewDM_->cbegin();
    auto it_Tau_rawMVAoldDM_ = Tau_rawMVAoldDM_->cbegin();
    auto it_Tau_rawMVAoldDMdR03_ = Tau_rawMVAoldDMdR03_->cbegin();
    auto it_Tau_charge_ = Tau_charge_->cbegin();
    auto it_Tau_decayMode_ = Tau_decayMode_->cbegin();
    auto it_Tau_jetIdx_ = Tau_jetIdx_->cbegin();
    auto it_Tau_rawAntiEleCat_ = Tau_rawAntiEleCat_->cbegin();
    auto it_Tau_idAntiEle_ = Tau_idAntiEle_->cbegin();
    auto it_Tau_idAntiMu_ = Tau_idAntiMu_->cbegin();
    auto it_Tau_idDecayMode_ = Tau_idDecayMode_->cbegin();
    auto it_Tau_idDecayModeNewDMs_ = Tau_idDecayModeNewDMs_->cbegin();
    auto it_Tau_idMVAnewDM_ = Tau_idMVAnewDM_->cbegin();
    auto it_Tau_idMVAoldDM_ = Tau_idMVAoldDM_->cbegin();
    auto it_Tau_idMVAoldDMdR03_ = Tau_idMVAoldDMdR03_->cbegin();
    auto it_Tau_cleanmask_ = Tau_cleanmask_->cbegin();
    auto it_Tau_genPartIdx_ = Tau_genPartIdx_->cbegin();
    auto it_Tau_genPartFlav_ = Tau_genPartFlav_->cbegin();
    for(; it_Tau_chargedIso_ != Tau_chargedIso_->cend(); ){
      Ta obj;
      obj.setchargedIso(*it_Tau_chargedIso_);
      obj.setdxy(*it_Tau_dxy_);
      obj.setdz(*it_Tau_dz_);
      obj.setfootprintCorr(*it_Tau_footprintCorr_);
      obj.setleadTkDeltaEta(*it_Tau_leadTkDeltaEta_);
      obj.setleadTkDeltaPhi(*it_Tau_leadTkDeltaPhi_);
      obj.setleadTkPtOverTauPt(*it_Tau_leadTkPtOverTauPt_);
      obj.setneutralIso(*it_Tau_neutralIso_);
      obj.setphotonsOutsideSignalCone(*it_Tau_photonsOutsideSignalCone_);
      obj.setpuCorr(*it_Tau_puCorr_);
      obj.setrawAntiEle(*it_Tau_rawAntiEle_);
      obj.setrawIso(*it_Tau_rawIso_);
      obj.setrawMVAnewDM(*it_Tau_rawMVAnewDM_);
      obj.setrawMVAoldDM(*it_Tau_rawMVAoldDM_);
      obj.setrawMVAoldDMdR03(*it_Tau_rawMVAoldDMdR03_);
      obj.setcharge(*it_Tau_charge_);
      obj.setdecayMode(*it_Tau_decayMode_);
      obj.setjetIdx(*it_Tau_jetIdx_);
      obj.setrawAntiEleCat(*it_Tau_rawAntiEleCat_);
      obj.setidAntiEle(*it_Tau_idAntiEle_);
      obj.setidAntiMu(*it_Tau_idAntiMu_);
      obj.setidDecayMode(*it_Tau_idDecayMode_);
      obj.setidDecayModeNewDMs(*it_Tau_idDecayModeNewDMs_);
      obj.setidMVAnewDM(*it_Tau_idMVAnewDM_);
      obj.setidMVAoldDM(*it_Tau_idMVAoldDM_);
      obj.setidMVAoldDMdR03(*it_Tau_idMVAoldDMdR03_);
      obj.setcleanmask(*it_Tau_cleanmask_);
      obj.setgenPartIdx(*it_Tau_genPartIdx_);
      obj.setgenPartFlav(*it_Tau_genPartFlav_);
      obj.setLotentzVector(*it_Tau_pt_, *it_Tau_eta_, *it_Tau_phi_, *it_Tau_mass_);
      Tau_.push_back( obj );
      ++it_Tau_chargedIso_;
      ++it_Tau_dxy_;
      ++it_Tau_dz_;
      ++it_Tau_eta_;
      ++it_Tau_footprintCorr_;
      ++it_Tau_leadTkDeltaEta_;
      ++it_Tau_leadTkDeltaPhi_;
      ++it_Tau_leadTkPtOverTauPt_;
      ++it_Tau_mass_;
      ++it_Tau_neutralIso_;
      ++it_Tau_phi_;
      ++it_Tau_photonsOutsideSignalCone_;
      ++it_Tau_pt_;
      ++it_Tau_puCorr_;
      ++it_Tau_rawAntiEle_;
      ++it_Tau_rawIso_;
      ++it_Tau_rawMVAnewDM_;
      ++it_Tau_rawMVAoldDM_;
      ++it_Tau_rawMVAoldDMdR03_;
      ++it_Tau_charge_;
      ++it_Tau_decayMode_;
      ++it_Tau_jetIdx_;
      ++it_Tau_rawAntiEleCat_;
      ++it_Tau_idAntiEle_;
      ++it_Tau_idAntiMu_;
      ++it_Tau_idDecayMode_;
      ++it_Tau_idDecayModeNewDMs_;
      ++it_Tau_idMVAnewDM_;
      ++it_Tau_idMVAoldDM_;
      ++it_Tau_idMVAoldDMdR03_;
      ++it_Tau_cleanmask_;
      ++it_Tau_genPartIdx_;
      ++it_Tau_genPartFlav_;
    }
    return Tau_;
  }
  
  const Puppimet PuppiMET(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadPuppimet();
  
    Puppimet obj;
    obj.setsumEt(PuppiMET_sumEt_);
  
    return obj;
  }
  
  const vector<Muo>& Muon(){
    if(Muon_.size() > 0) return Muon_;
    loadMuos();
  	Muon_.reserve(Muon_dxy_->size());
    auto it_Muon_dxy_ = Muon_dxy_->cbegin();
    auto it_Muon_dxyErr_ = Muon_dxyErr_->cbegin();
    auto it_Muon_dz_ = Muon_dz_->cbegin();
    auto it_Muon_dzErr_ = Muon_dzErr_->cbegin();
    auto it_Muon_eta_ = Muon_eta_->cbegin();
    auto it_Muon_ip3d_ = Muon_ip3d_->cbegin();
    auto it_Muon_mass_ = Muon_mass_->cbegin();
    auto it_Muon_miniPFRelIso_all_ = Muon_miniPFRelIso_all_->cbegin();
    auto it_Muon_miniPFRelIso_chg_ = Muon_miniPFRelIso_chg_->cbegin();
    auto it_Muon_pfRelIso03_all_ = Muon_pfRelIso03_all_->cbegin();
    auto it_Muon_pfRelIso03_chg_ = Muon_pfRelIso03_chg_->cbegin();
    auto it_Muon_pfRelIso04_all_ = Muon_pfRelIso04_all_->cbegin();
    auto it_Muon_phi_ = Muon_phi_->cbegin();
    auto it_Muon_pt_ = Muon_pt_->cbegin();
    auto it_Muon_ptErr_ = Muon_ptErr_->cbegin();
    auto it_Muon_segmentComp_ = Muon_segmentComp_->cbegin();
    auto it_Muon_sip3d_ = Muon_sip3d_->cbegin();
    auto it_Muon_mvaTTH_ = Muon_mvaTTH_->cbegin();
    auto it_Muon_charge_ = Muon_charge_->cbegin();
    auto it_Muon_jetIdx_ = Muon_jetIdx_->cbegin();
    auto it_Muon_nStations_ = Muon_nStations_->cbegin();
    auto it_Muon_nTrackerLayers_ = Muon_nTrackerLayers_->cbegin();
    auto it_Muon_pdgId_ = Muon_pdgId_->cbegin();
    auto it_Muon_tightCharge_ = Muon_tightCharge_->cbegin();
    auto it_Muon_highPtId_ = Muon_highPtId_->cbegin();
    auto it_Muon_isPFcand_ = Muon_isPFcand_->cbegin();
    auto it_Muon_mediumId_ = Muon_mediumId_->cbegin();
    auto it_Muon_softId_ = Muon_softId_->cbegin();
    auto it_Muon_tightId_ = Muon_tightId_->cbegin();
    auto it_Muon_genPartIdx_ = Muon_genPartIdx_->cbegin();
    auto it_Muon_genPartFlav_ = Muon_genPartFlav_->cbegin();
    auto it_Muon_cleanmask_ = Muon_cleanmask_->cbegin();
    for(; it_Muon_dxy_ != Muon_dxy_->cend(); ){
      Muo obj;
      obj.setdxy(*it_Muon_dxy_);
      obj.setdxyErr(*it_Muon_dxyErr_);
      obj.setdz(*it_Muon_dz_);
      obj.setdzErr(*it_Muon_dzErr_);
      obj.setip3d(*it_Muon_ip3d_);
      obj.setminiPFRelIso(*it_Muon_miniPFRelIso_all_);
      obj.setminiPFRelIso(*it_Muon_miniPFRelIso_chg_);
      obj.setpfRelIso03(*it_Muon_pfRelIso03_all_);
      obj.setpfRelIso03(*it_Muon_pfRelIso03_chg_);
      obj.setpfRelIso04(*it_Muon_pfRelIso04_all_);
      obj.setptErr(*it_Muon_ptErr_);
      obj.setsegmentComp(*it_Muon_segmentComp_);
      obj.setsip3d(*it_Muon_sip3d_);
      obj.setmvaTTH(*it_Muon_mvaTTH_);
      obj.setcharge(*it_Muon_charge_);
      obj.setjetIdx(*it_Muon_jetIdx_);
      obj.setnStations(*it_Muon_nStations_);
      obj.setnTrackerLayers(*it_Muon_nTrackerLayers_);
      obj.setpdgId(*it_Muon_pdgId_);
      obj.settightCharge(*it_Muon_tightCharge_);
      obj.sethighPtId(*it_Muon_highPtId_);
      obj.setisPFcand(*it_Muon_isPFcand_);
      obj.setmediumId(*it_Muon_mediumId_);
      obj.setsoftId(*it_Muon_softId_);
      obj.settightId(*it_Muon_tightId_);
      obj.setgenPartIdx(*it_Muon_genPartIdx_);
      obj.setgenPartFlav(*it_Muon_genPartFlav_);
      obj.setcleanmask(*it_Muon_cleanmask_);
      obj.setLotentzVector(*it_Muon_pt_, *it_Muon_eta_, *it_Muon_phi_, *it_Muon_mass_);
      Muon_.push_back( obj );
      ++it_Muon_dxy_;
      ++it_Muon_dxyErr_;
      ++it_Muon_dz_;
      ++it_Muon_dzErr_;
      ++it_Muon_eta_;
      ++it_Muon_ip3d_;
      ++it_Muon_mass_;
      ++it_Muon_miniPFRelIso_all_;
      ++it_Muon_miniPFRelIso_chg_;
      ++it_Muon_pfRelIso03_all_;
      ++it_Muon_pfRelIso03_chg_;
      ++it_Muon_pfRelIso04_all_;
      ++it_Muon_phi_;
      ++it_Muon_pt_;
      ++it_Muon_ptErr_;
      ++it_Muon_segmentComp_;
      ++it_Muon_sip3d_;
      ++it_Muon_mvaTTH_;
      ++it_Muon_charge_;
      ++it_Muon_jetIdx_;
      ++it_Muon_nStations_;
      ++it_Muon_nTrackerLayers_;
      ++it_Muon_pdgId_;
      ++it_Muon_tightCharge_;
      ++it_Muon_highPtId_;
      ++it_Muon_isPFcand_;
      ++it_Muon_mediumId_;
      ++it_Muon_softId_;
      ++it_Muon_tightId_;
      ++it_Muon_genPartIdx_;
      ++it_Muon_genPartFlav_;
      ++it_Muon_cleanmask_;
    }
    return Muon_;
  }
  
  const vector<Otherp>& OtherPV(){
    if(OtherPV_.size() > 0) return OtherPV_;
    loadOtherps();
  	OtherPV_.reserve(OtherPV_z_->size());
    auto it_OtherPV_z_ = OtherPV_z_->cbegin();
    for(; it_OtherPV_z_ != OtherPV_z_->cend(); ){
      Otherp obj;
      obj.setz(*it_OtherPV_z_);
      
      OtherPV_.push_back( obj );
      ++it_OtherPV_z_;
    }
    return OtherPV_;
  }
  
  const Hlt HLT(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadHlt();
  
    Hlt obj;
    obj.setHLTriggerFirstPath(HLTriggerFirstPath);
    obj.setAK8PFJet360(HLT_AK8PFJet360_TrimMass30_);
    obj.setAK8PFJet400(HLT_AK8PFJet400_TrimMass30_);
    obj.setAK8PFHT750(HLT_AK8PFHT750_TrimMass50_);
    obj.setAK8PFHT800(HLT_AK8PFHT800_TrimMass50_);
    obj.setAK8DiPFJet300(HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_);
    obj.setAK8DiPFJet280(HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_);
    obj.setAK8DiPFJet300(HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_);
    obj.setAK8DiPFJet300(HLT_AK8DiPFJet300_200_TrimMass30_);
    obj.setAK8PFHT700(HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_);
    obj.setAK8PFHT650(HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_);
    obj.setAK8PFHT600(HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_);
    obj.setAK8DiPFJet280(HLT_AK8DiPFJet280_200_TrimMass30_);
    obj.setAK8DiPFJet250(HLT_AK8DiPFJet250_200_TrimMass30_);
    obj.setAK8DiPFJet280(HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_);
    obj.setAK8DiPFJet250(HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_);
    obj.setCaloJet260(HLT_CaloJet260_);
    obj.setCaloJet500(HLT_CaloJet500_NoJetID_);
    obj.setDimuon13(HLT_Dimuon13_PsiPrime_);
    obj.setDimuon13(HLT_Dimuon13_Upsilon_);
    obj.setDimuon20(HLT_Dimuon20_Jpsi_);
    obj.setDoubleEle24(HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_);
    obj.setDoubleEle25(HLT_DoubleEle25_CaloIdL_GsfTrkIdVL_);
    obj.setDoubleEle33(HLT_DoubleEle33_CaloIdL_);
    obj.setDoubleEle33(HLT_DoubleEle33_CaloIdL_MW_);
    obj.setDoubleEle33(HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_);
    obj.setDoubleEle33(HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_);
    obj.setDoubleMediumCombinedIsoPFTau35(HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_);
    obj.setDoubleTightCombinedIsoPFTau35(HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg_);
    obj.setDoubleMediumCombinedIsoPFTau40(HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_);
    obj.setDoubleTightCombinedIsoPFTau40(HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg_);
    obj.setDoubleMediumCombinedIsoPFTau40(HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_);
    obj.setDoubleTightCombinedIsoPFTau40(HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_);
    obj.setDoubleMediumIsoPFTau35(HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_);
    obj.setDoubleMediumIsoPFTau40(HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_);
    obj.setDoubleMediumIsoPFTau40(HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_);
    obj.setDoubleEle37(HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_);
    obj.setDoubleMu33NoFiltersNoVtx(HLT_DoubleMu33NoFiltersNoVtx_);
    obj.setDoubleMu38NoFiltersNoVtx(HLT_DoubleMu38NoFiltersNoVtx_);
    obj.setDoubleMu23NoFiltersNoVtxDisplaced(HLT_DoubleMu23NoFiltersNoVtxDisplaced_);
    obj.setDoubleMu28NoFiltersNoVtxDisplaced(HLT_DoubleMu28NoFiltersNoVtxDisplaced_);
    obj.setDoubleMu0(HLT_DoubleMu0_);
    obj.setDoubleMu4(HLT_DoubleMu4_3_Bs_);
    obj.setDoubleMu4(HLT_DoubleMu4_3_Jpsi_Displaced_);
    obj.setDoubleMu4(HLT_DoubleMu4_JpsiTrk_Displaced_);
    obj.setDoubleMu4(HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_);
    obj.setDoubleMu3(HLT_DoubleMu3_Trk_Tau3mu_);
    obj.setDoubleMu4(HLT_DoubleMu4_PsiPrimeTrk_Displaced_);
    obj.setMu7p5(HLT_Mu7p5_L2Mu2_Jpsi_);
    obj.setMu7p5(HLT_Mu7p5_L2Mu2_Upsilon_);
    obj.setMu7p5(HLT_Mu7p5_Track2_Jpsi_);
    obj.setMu7p5(HLT_Mu7p5_Track3p5_Jpsi_);
    obj.setMu7p5(HLT_Mu7p5_Track7_Jpsi_);
    obj.setMu7p5(HLT_Mu7p5_Track2_Upsilon_);
    obj.setMu7p5(HLT_Mu7p5_Track3p5_Upsilon_);
    obj.setMu7p5(HLT_Mu7p5_Track7_Upsilon_);
    obj.setDimuon0er16(HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_);
    obj.setDimuon0er16(HLT_Dimuon0er16_Jpsi_NoVertexing_);
    obj.setDimuon6(HLT_Dimuon6_Jpsi_NoVertexing_);
    obj.setPhoton150(HLT_Photon150_);
    obj.setPhoton90(HLT_Photon90_CaloIdL_HT300_);
    obj.setHT250(HLT_HT250_CaloMET70_);
    obj.setDoublePhoton60(HLT_DoublePhoton60_);
    obj.setDoublePhoton85(HLT_DoublePhoton85_);
    obj.setEle17(HLT_Ele17_Ele8_Gsf_);
    obj.setEle20(HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28_);
    obj.setEle22(HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29_);
    obj.setEle22(HLT_Ele22_eta2p1_WPLoose_Gsf_);
    obj.setEle22(HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
    obj.setEle23(HLT_Ele23_WPLoose_Gsf_);
    obj.setEle23(HLT_Ele23_WPLoose_Gsf_WHbbBoost_);
    obj.setEle24(HLT_Ele24_eta2p1_WPLoose_Gsf_);
    obj.setEle24(HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_);
    obj.setEle24(HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
    obj.setEle24(HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_);
    obj.setEle25(HLT_Ele25_WPTight_Gsf_);
    obj.setEle25(HLT_Ele25_eta2p1_WPLoose_Gsf_);
    obj.setEle25(HLT_Ele25_eta2p1_WPTight_Gsf_);
    obj.setEle27(HLT_Ele27_WPLoose_Gsf_);
    obj.setEle27(HLT_Ele27_WPLoose_Gsf_WHbbBoost_);
    obj.setEle27(HLT_Ele27_WPTight_Gsf_);
    obj.setEle27(HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_);
    obj.setEle27(HLT_Ele27_eta2p1_WPLoose_Gsf_);
    obj.setEle27(HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
    obj.setEle27(HLT_Ele27_eta2p1_WPTight_Gsf_);
    obj.setEle30(HLT_Ele30_WPTight_Gsf_);
    obj.setEle30(HLT_Ele30_eta2p1_WPLoose_Gsf_);
    obj.setEle30(HLT_Ele30_eta2p1_WPTight_Gsf_);
    obj.setEle32(HLT_Ele32_WPTight_Gsf_);
    obj.setEle32(HLT_Ele32_eta2p1_WPLoose_Gsf_);
    obj.setEle32(HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
    obj.setEle32(HLT_Ele32_eta2p1_WPTight_Gsf_);
    obj.setEle35(HLT_Ele35_WPLoose_Gsf_);
    obj.setEle35(HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_);
    obj.setEle36(HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_);
    obj.setEle45(HLT_Ele45_WPLoose_Gsf_);
    obj.setEle45(HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded_);
    obj.setEle45(HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_);
    obj.setEle105(HLT_Ele105_CaloIdVT_GsfTrkIdT_);
    obj.setEle30WP60(HLT_Ele30WP60_SC4_Mass55_);
    obj.setEle30WP60(HLT_Ele30WP60_Ele8_Mass55_);
    obj.setHT200(HLT_HT200_);
    obj.setHT275(HLT_HT275_);
    obj.setHT325(HLT_HT325_);
    obj.setHT425(HLT_HT425_);
    obj.setHT575(HLT_HT575_);
    obj.setHT410to430(HLT_HT410to430_);
    obj.setHT430to450(HLT_HT430to450_);
    obj.setHT450to470(HLT_HT450to470_);
    obj.setHT470to500(HLT_HT470to500_);
    obj.setHT500to550(HLT_HT500to550_);
    obj.setHT550to650(HLT_HT550to650_);
    obj.setHT650(HLT_HT650_);
    obj.setMu16(HLT_Mu16_eta2p1_MET30_);
    obj.setIsoMu16(HLT_IsoMu16_eta2p1_MET30_);
    obj.setIsoMu16(HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1_);
    obj.setIsoMu17(HLT_IsoMu17_eta2p1_);
    obj.setIsoMu17(HLT_IsoMu17_eta2p1_LooseIsoPFTau20_);
    obj.setIsoMu17(HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_);
    obj.setDoubleIsoMu17(HLT_DoubleIsoMu17_eta2p1_);
    obj.setDoubleIsoMu17(HLT_DoubleIsoMu17_eta2p1_noDzCut_);
    obj.setIsoMu18(HLT_IsoMu18_);
    obj.setIsoMu19(HLT_IsoMu19_eta2p1_LooseIsoPFTau20_);
    obj.setIsoMu19(HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_);
    obj.setIsoMu19(HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_);
    obj.setIsoMu19(HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20_);
    obj.setIsoMu19(HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_);
    obj.setIsoMu19(HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_);
    obj.setIsoMu21(HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_);
    obj.setIsoMu21(HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg_);
    obj.setIsoMu20(HLT_IsoMu20_);
    obj.setIsoMu21(HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_);
    obj.setIsoMu21(HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1_);
    obj.setIsoMu21(HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_);
    obj.setIsoMu22(HLT_IsoMu22_);
    obj.setIsoMu22(HLT_IsoMu22_eta2p1_);
    obj.setIsoMu24(HLT_IsoMu24_);
    obj.setIsoMu27(HLT_IsoMu27_);
    obj.setIsoTkMu18(HLT_IsoTkMu18_);
    obj.setIsoTkMu20(HLT_IsoTkMu20_);
    obj.setIsoTkMu22(HLT_IsoTkMu22_);
    obj.setIsoTkMu22(HLT_IsoTkMu22_eta2p1_);
    obj.setIsoTkMu24(HLT_IsoTkMu24_);
    obj.setIsoTkMu27(HLT_IsoTkMu27_);
    obj.setJetE30(HLT_JetE30_NoBPTX3BX_);
    obj.setJetE30(HLT_JetE30_NoBPTX_);
    obj.setJetE50(HLT_JetE50_NoBPTX3BX_);
    obj.setJetE70(HLT_JetE70_NoBPTX3BX_);
    obj.setL1SingleMu18(HLT_L1SingleMu18_);
    obj.setL2Mu10(HLT_L2Mu10_);
    obj.setL1SingleMuOpen(HLT_L1SingleMuOpen_);
    obj.setL1SingleMuOpen(HLT_L1SingleMuOpen_DT_);
    obj.setL2DoubleMu23(HLT_L2DoubleMu23_NoVertex_);
    obj.setL2DoubleMu28(HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_);
    obj.setL2DoubleMu38(HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10_);
    obj.setL2Mu10(HLT_L2Mu10_NoVertex_NoBPTX3BX_);
    obj.setL2Mu10(HLT_L2Mu10_NoVertex_NoBPTX_);
    obj.setL2Mu45(HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_);
    obj.setL2Mu40(HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_);
    obj.setLooseIsoPFTau50(HLT_LooseIsoPFTau50_Trk30_eta2p1_);
    obj.setLooseIsoPFTau50(HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_);
    obj.setLooseIsoPFTau50(HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_);
    obj.setLooseIsoPFTau50(HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_);
    obj.setLooseIsoPFTau50(HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_);
    obj.setPFTau120(HLT_PFTau120_eta2p1_);
    obj.setPFTau140(HLT_PFTau140_eta2p1_);
    obj.setVLooseIsoPFTau120(HLT_VLooseIsoPFTau120_Trk50_eta2p1_);
    obj.setVLooseIsoPFTau140(HLT_VLooseIsoPFTau140_Trk50_eta2p1_);
    obj.setMu17(HLT_Mu17_Mu8_);
    obj.setMu17(HLT_Mu17_Mu8_DZ_);
    obj.setMu17(HLT_Mu17_Mu8_SameSign_);
    obj.setMu17(HLT_Mu17_Mu8_SameSign_DZ_);
    obj.setMu20(HLT_Mu20_Mu10_);
    obj.setMu20(HLT_Mu20_Mu10_DZ_);
    obj.setMu20(HLT_Mu20_Mu10_SameSign_);
    obj.setMu20(HLT_Mu20_Mu10_SameSign_DZ_);
    obj.setMu17(HLT_Mu17_TkMu8_DZ_);
    obj.setMu17(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_);
    obj.setMu17(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_);
    obj.setMu17(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_);
    obj.setMu17(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_);
    obj.setMu25(HLT_Mu25_TkMu0_dEta18_Onia_);
    obj.setMu27(HLT_Mu27_TkMu8_);
    obj.setMu30(HLT_Mu30_TkMu11_);
    obj.setMu30(HLT_Mu30_eta2p1_PFJet150_PFJet50_);
    obj.setMu40(HLT_Mu40_TkMu11_);
    obj.setMu40(HLT_Mu40_eta2p1_PFJet200_PFJet50_);
    obj.setMu20(HLT_Mu20_);
    obj.setTkMu17(HLT_TkMu17_);
    obj.setTkMu17(HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_);
    obj.setTkMu17(HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_);
    obj.setTkMu20(HLT_TkMu20_);
    obj.setMu24(HLT_Mu24_eta2p1_);
    obj.setTkMu24(HLT_TkMu24_eta2p1_);
    obj.setMu27(HLT_Mu27_);
    obj.setTkMu27(HLT_TkMu27_);
    obj.setMu45(HLT_Mu45_eta2p1_);
    obj.setMu50(HLT_Mu50_);
    obj.setTkMu50(HLT_TkMu50_);
    obj.setMu38NoFiltersNoVtx(HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_);
    obj.setMu42NoFiltersNoVtx(HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL_);
    obj.setMu28NoFiltersNoVtxDisplaced(HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL_);
    obj.setMu33NoFiltersNoVtxDisplaced(HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL_);
    obj.setMu23NoFiltersNoVtx(HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL_);
    obj.setDoubleMu18NoFiltersNoVtx(HLT_DoubleMu18NoFiltersNoVtx_);
    obj.setMu33NoFiltersNoVtxDisplaced(HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight_);
    obj.setMu33NoFiltersNoVtxDisplaced(HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose_);
    obj.setMu28NoFiltersNoVtx(HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose_);
    obj.setMu38NoFiltersNoVtxDisplaced(HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight_);
    obj.setMu38NoFiltersNoVtxDisplaced(HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose_);
    obj.setMu38NoFiltersNoVtx(HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose_);
    obj.setMu28NoFiltersNoVtx(HLT_Mu28NoFiltersNoVtx_CentralCaloJet40_);
    obj.setPFHT300(HLT_PFHT300_PFMET100_);
    obj.setPFHT300(HLT_PFHT300_PFMET110_);
    obj.setPFHT550(HLT_PFHT550_4JetPt50_);
    obj.setPFHT650(HLT_PFHT650_4JetPt50_);
    obj.setPFHT750(HLT_PFHT750_4JetPt50_);
    obj.setPFHT750(HLT_PFHT750_4JetPt70_);
    obj.setPFHT750(HLT_PFHT750_4JetPt80_);
    obj.setPFHT800(HLT_PFHT800_4JetPt50_);
    obj.setPFHT850(HLT_PFHT850_4JetPt50_);
    obj.setPFJet15(HLT_PFJet15_NoCaloMatched_);
    obj.setPFJet25(HLT_PFJet25_NoCaloMatched_);
    obj.setDiPFJet15(HLT_DiPFJet15_NoCaloMatched_);
    obj.setDiPFJet25(HLT_DiPFJet25_NoCaloMatched_);
    obj.setDiPFJet15(HLT_DiPFJet15_FBEta3_NoCaloMatched_);
    obj.setDiPFJet25(HLT_DiPFJet25_FBEta3_NoCaloMatched_);
    obj.setDiPFJetAve15(HLT_DiPFJetAve15_HFJEC_);
    obj.setDiPFJetAve25(HLT_DiPFJetAve25_HFJEC_);
    obj.setDiPFJetAve35(HLT_DiPFJetAve35_HFJEC_);
    obj.setAK8PFJet40(HLT_AK8PFJet40_);
    obj.setAK8PFJet60(HLT_AK8PFJet60_);
    obj.setAK8PFJet80(HLT_AK8PFJet80_);
    obj.setAK8PFJet140(HLT_AK8PFJet140_);
    obj.setAK8PFJet200(HLT_AK8PFJet200_);
    obj.setAK8PFJet260(HLT_AK8PFJet260_);
    obj.setAK8PFJet320(HLT_AK8PFJet320_);
    obj.setAK8PFJet400(HLT_AK8PFJet400_);
    obj.setAK8PFJet450(HLT_AK8PFJet450_);
    obj.setAK8PFJet500(HLT_AK8PFJet500_);
    obj.setPFJet40(HLT_PFJet40_);
    obj.setPFJet60(HLT_PFJet60_);
    obj.setPFJet80(HLT_PFJet80_);
    obj.setPFJet140(HLT_PFJet140_);
    obj.setPFJet200(HLT_PFJet200_);
    obj.setPFJet260(HLT_PFJet260_);
    obj.setPFJet320(HLT_PFJet320_);
    obj.setPFJet400(HLT_PFJet400_);
    obj.setPFJet450(HLT_PFJet450_);
    obj.setPFJet500(HLT_PFJet500_);
    obj.setDiPFJetAve40(HLT_DiPFJetAve40_);
    obj.setDiPFJetAve60(HLT_DiPFJetAve60_);
    obj.setDiPFJetAve80(HLT_DiPFJetAve80_);
    obj.setDiPFJetAve140(HLT_DiPFJetAve140_);
    obj.setDiPFJetAve200(HLT_DiPFJetAve200_);
    obj.setDiPFJetAve260(HLT_DiPFJetAve260_);
    obj.setDiPFJetAve320(HLT_DiPFJetAve320_);
    obj.setDiPFJetAve400(HLT_DiPFJetAve400_);
    obj.setDiPFJetAve500(HLT_DiPFJetAve500_);
    obj.setDiPFJetAve60(HLT_DiPFJetAve60_HFJEC_);
    obj.setDiPFJetAve80(HLT_DiPFJetAve80_HFJEC_);
    obj.setDiPFJetAve100(HLT_DiPFJetAve100_HFJEC_);
    obj.setDiPFJetAve160(HLT_DiPFJetAve160_HFJEC_);
    obj.setDiPFJetAve220(HLT_DiPFJetAve220_HFJEC_);
    obj.setDiPFJetAve300(HLT_DiPFJetAve300_HFJEC_);
    obj.setDiPFJet40(HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140_);
    obj.setDiPFJet40(HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80_);
    obj.setDiCentralPFJet170(HLT_DiCentralPFJet170_);
    obj.setSingleCentralPFJet170(HLT_SingleCentralPFJet170_CFMax0p1_);
    obj.setDiCentralPFJet170(HLT_DiCentralPFJet170_CFMax0p1_);
    obj.setDiCentralPFJet220(HLT_DiCentralPFJet220_CFMax0p3_);
    obj.setDiCentralPFJet330(HLT_DiCentralPFJet330_CFMax0p5_);
    obj.setDiCentralPFJet430(HLT_DiCentralPFJet430_);
    obj.setPFHT125(HLT_PFHT125_);
    obj.setPFHT200(HLT_PFHT200_);
    obj.setPFHT250(HLT_PFHT250_);
    obj.setPFHT300(HLT_PFHT300_);
    obj.setPFHT350(HLT_PFHT350_);
    obj.setPFHT400(HLT_PFHT400_);
    obj.setPFHT475(HLT_PFHT475_);
    obj.setPFHT600(HLT_PFHT600_);
    obj.setPFHT650(HLT_PFHT650_);
    obj.setPFHT800(HLT_PFHT800_);
    obj.setPFHT900(HLT_PFHT900_);
    obj.setPFHT200(HLT_PFHT200_PFAlphaT0p51_);
    obj.setPFHT200(HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57_);
    obj.setPFHT200(HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_);
    obj.setPFHT250(HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55_);
    obj.setPFHT250(HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_);
    obj.setPFHT300(HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53_);
    obj.setPFHT300(HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_);
    obj.setPFHT350(HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52_);
    obj.setPFHT350(HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53_);
    obj.setPFHT400(HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51_);
    obj.setPFHT400(HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52_);
    obj.setMET60(HLT_MET60_IsoTrk35_Loose_);
    obj.setMET75(HLT_MET75_IsoTrk50_);
    obj.setMET90(HLT_MET90_IsoTrk50_);
    obj.setPFMET120(HLT_PFMET120_BTagCSV_p067_);
    obj.setPFMET120(HLT_PFMET120_Mu5_);
    obj.setPFMET170(HLT_PFMET170_NotCleaned_);
    obj.setPFMET170(HLT_PFMET170_NoiseCleaned_);
    obj.setPFMET170(HLT_PFMET170_HBHECleaned_);
    obj.setPFMET170(HLT_PFMET170_JetIdCleaned_);
    obj.setPFMET170(HLT_PFMET170_BeamHaloCleaned_);
    obj.setPFMET170(HLT_PFMET170_HBHE_BeamHaloCleaned_);
    obj.setPFMETTypeOne190(HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned_);
    obj.setPFMET90(HLT_PFMET90_PFMHT90_IDTight_);
    obj.setPFMET100(HLT_PFMET100_PFMHT100_IDTight_);
    obj.setPFMET100(HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned_);
    obj.setPFMET110(HLT_PFMET110_PFMHT110_IDTight_);
    obj.setPFMET120(HLT_PFMET120_PFMHT120_IDTight_);
    obj.setCaloMHTNoPU90(HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_);
    obj.setCaloMHTNoPU90(HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_);
    obj.setQuadPFJet(HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200_);
    obj.setQuadPFJet(HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460_);
    obj.setQuadPFJet(HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_);
    obj.setQuadPFJet(HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_);
    obj.setQuadPFJet(HLT_QuadPFJet_VBF_);
    obj.setL1(HLT_L1_TripleJet_VBF_);
    obj.setQuadJet45(HLT_QuadJet45_TripleBTagCSV_p087_);
    obj.setQuadJet45(HLT_QuadJet45_DoubleBTagCSV_p087_);
    obj.setDoubleJet90(HLT_DoubleJet90_Double30_TripleBTagCSV_p087_);
    obj.setDoubleJet90(HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_);
    obj.setDoubleJetsC100(HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_);
    obj.setDoubleJetsC100(HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_);
    obj.setDoubleJetsC112(HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172_);
    obj.setDoubleJetsC112(HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_);
    obj.setDoubleJetsC100(HLT_DoubleJetsC100_SingleBTagCSV_p026_);
    obj.setDoubleJetsC100(HLT_DoubleJetsC100_SingleBTagCSV_p014_);
    obj.setDoubleJetsC100(HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350_);
    obj.setDoubleJetsC100(HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350_);
    obj.setPhoton135(HLT_Photon135_PFMET100_);
    obj.setPhoton20(HLT_Photon20_CaloIdVL_IsoL_);
    obj.setPhoton22(HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
    obj.setPhoton22(HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_);
    obj.setPhoton250(HLT_Photon250_NoHE_);
    obj.setPhoton300(HLT_Photon300_NoHE_);
    obj.setPhoton26(HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60_);
    obj.setPhoton36(HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15_);
    obj.setPhoton36(HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
    obj.setPhoton36(HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_);
    obj.setPhoton50(HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
    obj.setPhoton50(HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_);
    obj.setPhoton75(HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
    obj.setPhoton75(HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_);
    obj.setPhoton90(HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
    obj.setPhoton90(HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_);
    obj.setPhoton120(HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_);
    obj.setPhoton120(HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_);
    obj.setMu8(HLT_Mu8_TrkIsoVVL_);
    obj.setMu17(HLT_Mu17_TrkIsoVVL_);
    obj.setEle8(HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_);
    obj.setEle12(HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_);
    obj.setEle17(HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_);
    obj.setEle23(HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_);
    obj.setBTagMu(HLT_BTagMu_DiJet20_Mu5_);
    obj.setBTagMu(HLT_BTagMu_DiJet40_Mu5_);
    obj.setBTagMu(HLT_BTagMu_DiJet70_Mu5_);
    obj.setBTagMu(HLT_BTagMu_DiJet110_Mu5_);
    obj.setBTagMu(HLT_BTagMu_DiJet170_Mu5_);
    obj.setBTagMu(HLT_BTagMu_Jet300_Mu5_);
    obj.setBTagMu(HLT_BTagMu_AK8Jet300_Mu5_);
    obj.setEle23(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setEle23(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_);
    obj.setEle17(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setEle16(HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_);
    obj.setMu8(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu8(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu8(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setMu12(HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu12(HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setMu17(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu23(HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu23(HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setMu23(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu23(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setMu30(HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_);
    obj.setMu33(HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_);
    obj.setMu37(HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL_);
    obj.setMu27(HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL_);
    obj.setMu8(HLT_Mu8_DiEle12_CaloIdL_TrackIdL_);
    obj.setMu12(HLT_Mu12_Photon25_CaloIdL_);
    obj.setMu12(HLT_Mu12_Photon25_CaloIdL_L1ISO_);
    obj.setMu12(HLT_Mu12_Photon25_CaloIdL_L1OR_);
    obj.setMu17(HLT_Mu17_Photon22_CaloIdL_L1ISO_);
    obj.setMu17(HLT_Mu17_Photon30_CaloIdL_L1ISO_);
    obj.setMu17(HLT_Mu17_Photon35_CaloIdL_L1ISO_);
    obj.setDiMu9(HLT_DiMu9_Ele9_CaloIdL_TrackIdL_);
    obj.setTripleMu(HLT_TripleMu_5_3_3_);
    obj.setTripleMu(HLT_TripleMu_12_10_5_);
    obj.setMu3er(HLT_Mu3er_PFHT140_PFMET125_);
    obj.setMu6(HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067_);
    obj.setMu6(HLT_Mu6_PFHT200_PFMET100_);
    obj.setMu14er(HLT_Mu14er_PFMET100_);
    obj.setEle17(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_);
    obj.setEle23(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_);
    obj.setEle12(HLT_Ele12_CaloIdL_TrackIdL_IsoVL_);
    obj.setEle17(HLT_Ele17_CaloIdL_GsfTrkIdVL_);
    obj.setEle17(HLT_Ele17_CaloIdL_TrackIdL_IsoVL_);
    obj.setEle23(HLT_Ele23_CaloIdL_TrackIdL_IsoVL_);
    obj.setPFHT650(HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_);
    obj.setPFHT650(HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_);
    obj.setPhoton22(HLT_Photon22_);
    obj.setPhoton30(HLT_Photon30_);
    obj.setPhoton36(HLT_Photon36_);
    obj.setPhoton50(HLT_Photon50_);
    obj.setPhoton75(HLT_Photon75_);
    obj.setPhoton90(HLT_Photon90_);
    obj.setPhoton120(HLT_Photon120_);
    obj.setPhoton175(HLT_Photon175_);
    obj.setPhoton165(HLT_Photon165_HE10_);
    obj.setPhoton22(HLT_Photon22_R9Id90_HE10_IsoM_);
    obj.setPhoton30(HLT_Photon30_R9Id90_HE10_IsoM_);
    obj.setPhoton36(HLT_Photon36_R9Id90_HE10_IsoM_);
    obj.setPhoton50(HLT_Photon50_R9Id90_HE10_IsoM_);
    obj.setPhoton75(HLT_Photon75_R9Id90_HE10_IsoM_);
    obj.setPhoton90(HLT_Photon90_R9Id90_HE10_IsoM_);
    obj.setPhoton120(HLT_Photon120_R9Id90_HE10_IsoM_);
    obj.setPhoton165(HLT_Photon165_R9Id90_HE10_IsoM_);
    obj.setDiphoton30(HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_);
    obj.setDiphoton30(HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_);
    obj.setDiphoton30PV(HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_);
    obj.setDiphoton30(HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55_);
    obj.setDiphoton30EB(HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_);
    obj.setDimuon0(HLT_Dimuon0_Jpsi_Muon_);
    obj.setDimuon0(HLT_Dimuon0_Upsilon_Muon_);
    obj.setQuadMuon0(HLT_QuadMuon0_Dimuon0_Jpsi_);
    obj.setQuadMuon0(HLT_QuadMuon0_Dimuon0_Upsilon_);
    obj.setRsq0p25(HLT_Rsq0p25_Calo_);
    obj.setRsqMR240(HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo_);
    obj.setRsqMR240(HLT_RsqMR240_Rsq0p09_MR200_Calo_);
    obj.setRsq0p25(HLT_Rsq0p25_);
    obj.setRsq0p30(HLT_Rsq0p30_);
    obj.setRsqMR240(HLT_RsqMR240_Rsq0p09_MR200_);
    obj.setRsqMR240(HLT_RsqMR240_Rsq0p09_MR200_4jet_);
    obj.setRsqMR270(HLT_RsqMR270_Rsq0p09_MR200_);
    obj.setRsqMR270(HLT_RsqMR270_Rsq0p09_MR200_4jet_);
    obj.setRsq0p02(HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200_);
    obj.setRsq0p02(HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_);
    obj.setRsq0p02(HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_);
    obj.setRsq0p02(HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_);
    obj.setRsq0p02(HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200_);
    obj.setHT200(HLT_HT200_DisplacedDijet40_DisplacedTrack_);
    obj.setHT250(HLT_HT250_DisplacedDijet40_DisplacedTrack_);
    obj.setHT350(HLT_HT350_DisplacedDijet40_DisplacedTrack_);
    obj.setHT350(HLT_HT350_DisplacedDijet80_DisplacedTrack_);
    obj.setHT350(HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack_);
    obj.setHT350(HLT_HT350_DisplacedDijet40_Inclusive_);
    obj.setHT400(HLT_HT400_DisplacedDijet40_Inclusive_);
    obj.setHT500(HLT_HT500_DisplacedDijet40_Inclusive_);
    obj.setHT550(HLT_HT550_DisplacedDijet40_Inclusive_);
    obj.setHT550(HLT_HT550_DisplacedDijet80_Inclusive_);
    obj.setHT650(HLT_HT650_DisplacedDijet80_Inclusive_);
    obj.setHT750(HLT_HT750_DisplacedDijet80_Inclusive_);
    obj.setVBF(HLT_VBF_DisplacedJet40_DisplacedTrack_);
    obj.setVBF(HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5_);
    obj.setVBF(HLT_VBF_DisplacedJet40_TightID_DisplacedTrack_);
    obj.setVBF(HLT_VBF_DisplacedJet40_Hadronic_);
    obj.setVBF(HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack_);
    obj.setVBF(HLT_VBF_DisplacedJet40_TightID_Hadronic_);
    obj.setVBF(HLT_VBF_DisplacedJet40_VTightID_Hadronic_);
    obj.setVBF(HLT_VBF_DisplacedJet40_VVTightID_Hadronic_);
    obj.setVBF(HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_);
    obj.setVBF(HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_);
    obj.setPFMETNoMu90(HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_);
    obj.setPFMETNoMu100(HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_);
    obj.setPFMETNoMu110(HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_);
    obj.setPFMETNoMu120(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_);
    obj.setMonoCentralPFJet80(HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight_);
    obj.setMonoCentralPFJet80(HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight_);
    obj.setMonoCentralPFJet80(HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_);
    obj.setMonoCentralPFJet80(HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_);
    obj.setEle27(HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_);
    obj.setPhoton90(HLT_Photon90_CaloIdL_PFHT500_);
    obj.setDoubleMu8(HLT_DoubleMu8_Mass8_PFHT250_);
    obj.setMu8(HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250_);
    obj.setDoubleEle8(HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250_);
    obj.setDoubleMu8(HLT_DoubleMu8_Mass8_PFHT300_);
    obj.setMu8(HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_);
    obj.setDoubleEle8(HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_);
    obj.setMu10(HLT_Mu10_CentralPFJet30_BTagCSV_p13_);
    obj.setDoubleMu3(HLT_DoubleMu3_PFMET50_);
    obj.setEle10(HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_);
    obj.setEle15(HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_);
    obj.setEle15(HLT_Ele15_IsoVVVL_PFHT350_PFMET50_);
    obj.setEle15(HLT_Ele15_IsoVVVL_PFHT600_);
    obj.setEle15(HLT_Ele15_IsoVVVL_PFHT350_);
    obj.setEle15(HLT_Ele15_IsoVVVL_PFHT400_PFMET50_);
    obj.setEle15(HLT_Ele15_IsoVVVL_PFHT400_);
    obj.setEle50(HLT_Ele50_IsoVVVL_PFHT400_);
    obj.setMu8(HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_);
    obj.setMu10(HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_);
    obj.setMu15(HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_);
    obj.setMu15(HLT_Mu15_IsoVVVL_PFHT350_PFMET50_);
    obj.setMu15(HLT_Mu15_IsoVVVL_PFHT600_);
    obj.setMu15(HLT_Mu15_IsoVVVL_PFHT350_);
    obj.setMu15(HLT_Mu15_IsoVVVL_PFHT400_PFMET50_);
    obj.setMu15(HLT_Mu15_IsoVVVL_PFHT400_);
    obj.setMu50(HLT_Mu50_IsoVVVL_PFHT400_);
    obj.setDimuon16(HLT_Dimuon16_Jpsi_);
    obj.setDimuon10(HLT_Dimuon10_Jpsi_Barrel_);
    obj.setDimuon8(HLT_Dimuon8_PsiPrime_Barrel_);
    obj.setDimuon8(HLT_Dimuon8_Upsilon_Barrel_);
    obj.setDimuon0(HLT_Dimuon0_Phi_Barrel_);
    obj.setMu16(HLT_Mu16_TkMu0_dEta18_Onia_);
    obj.setMu16(HLT_Mu16_TkMu0_dEta18_Phi_);
    obj.setTrkMu15(HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_);
    obj.setTrkMu17(HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_);
    obj.setMu8(HLT_Mu8_);
    obj.setMu17(HLT_Mu17_);
    obj.setMu3(HLT_Mu3_PFJet40_);
    obj.setEle8(HLT_Ele8_CaloIdM_TrackIdM_PFJet30_);
    obj.setEle12(HLT_Ele12_CaloIdM_TrackIdM_PFJet30_);
    obj.setEle17(HLT_Ele17_CaloIdM_TrackIdM_PFJet30_);
    obj.setEle23(HLT_Ele23_CaloIdM_TrackIdM_PFJet30_);
    obj.setEle50(HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_);
    obj.setEle50(HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_);
    obj.setPFHT400(HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_);
    obj.setPFHT450(HLT_PFHT450_SixJet40_BTagCSV_p056_);
    obj.setPFHT400(HLT_PFHT400_SixJet30_);
    obj.setPFHT450(HLT_PFHT450_SixJet40_);
    obj.setEle115(HLT_Ele115_CaloIdVT_GsfTrkIdT_);
    obj.setMu55(HLT_Mu55_);
    obj.setPhoton42(HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_);
    obj.setPhoton90(HLT_Photon90_CaloIdL_PFHT600_);
    obj.setPixelTracks(HLT_PixelTracks_Multiplicity60ForEndOfFill_);
    obj.setPixelTracks(HLT_PixelTracks_Multiplicity85ForEndOfFill_);
    obj.setPixelTracks(HLT_PixelTracks_Multiplicity110ForEndOfFill_);
    obj.setPixelTracks(HLT_PixelTracks_Multiplicity135ForEndOfFill_);
    obj.setPixelTracks(HLT_PixelTracks_Multiplicity160ForEndOfFill_);
    obj.setFullTracks(HLT_FullTracks_Multiplicity80_);
    obj.setFullTracks(HLT_FullTracks_Multiplicity100_);
    obj.setFullTracks(HLT_FullTracks_Multiplicity130_);
    obj.setFullTracks(HLT_FullTracks_Multiplicity150_);
    obj.setECALHT800(HLT_ECALHT800_);
    obj.setDiSC30(HLT_DiSC30_18_EIso_AND_HE_Mass70_);
    obj.setPhoton125(HLT_Photon125_);
    obj.setMET100(HLT_MET100_);
    obj.setMET150(HLT_MET150_);
    obj.setMET200(HLT_MET200_);
    obj.setEle27(HLT_Ele27_HighEta_Ele20_Mass55_);
    obj.setL1FatEvents(HLT_L1FatEvents_);
    obj.setPhysics(HLT_Physics_);
    obj.setL1FatEvents(HLT_L1FatEvents_part0_);
    obj.setL1FatEvents(HLT_L1FatEvents_part1_);
    obj.setL1FatEvents(HLT_L1FatEvents_part2_);
    obj.setL1FatEvents(HLT_L1FatEvents_part3_);
    obj.setRandom(HLT_Random_);
    obj.setZeroBias(HLT_ZeroBias_);
    obj.setAK4CaloJet30(HLT_AK4CaloJet30_);
    obj.setAK4CaloJet40(HLT_AK4CaloJet40_);
    obj.setAK4CaloJet50(HLT_AK4CaloJet50_);
    obj.setAK4CaloJet80(HLT_AK4CaloJet80_);
    obj.setAK4CaloJet100(HLT_AK4CaloJet100_);
    obj.setAK4PFJet30(HLT_AK4PFJet30_);
    obj.setAK4PFJet50(HLT_AK4PFJet50_);
    obj.setAK4PFJet80(HLT_AK4PFJet80_);
    obj.setAK4PFJet100(HLT_AK4PFJet100_);
    obj.setHISinglePhoton10(HLT_HISinglePhoton10_);
    obj.setHISinglePhoton15(HLT_HISinglePhoton15_);
    obj.setHISinglePhoton20(HLT_HISinglePhoton20_);
    obj.setHISinglePhoton40(HLT_HISinglePhoton40_);
    obj.setHISinglePhoton60(HLT_HISinglePhoton60_);
    obj.setEcalCalibration(HLT_EcalCalibration_);
    obj.setHcalCalibration(HLT_HcalCalibration_);
    obj.setGlobalRunHPDNoise(HLT_GlobalRunHPDNoise_);
    obj.setL1BptxMinus(HLT_L1BptxMinus_);
    obj.setL1BptxPlus(HLT_L1BptxPlus_);
    obj.setL1NotBptxOR(HLT_L1NotBptxOR_);
    obj.setL1BeamGasMinus(HLT_L1BeamGasMinus_);
    obj.setL1BeamGasPlus(HLT_L1BeamGasPlus_);
    obj.setL1BptxXOR(HLT_L1BptxXOR_);
    obj.setL1MinimumBiasHF(HLT_L1MinimumBiasHF_OR_);
    obj.setL1MinimumBiasHF(HLT_L1MinimumBiasHF_AND_);
    obj.setHcalNZS(HLT_HcalNZS_);
    obj.setHcalPhiSym(HLT_HcalPhiSym_);
    obj.setHcalIsolatedbunch(HLT_HcalIsolatedbunch_);
    obj.setZeroBias(HLT_ZeroBias_FirstCollisionAfterAbortGap_);
    obj.setZeroBias(HLT_ZeroBias_FirstCollisionAfterAbortGap_copy_);
    obj.setZeroBias(HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS_);
    obj.setZeroBias(HLT_ZeroBias_IsolatedBunches_);
    obj.setZeroBias(HLT_ZeroBias_FirstCollisionInTrain_);
    obj.setZeroBias(HLT_ZeroBias_FirstBXAfterTrain_);
    obj.setPhoton500(HLT_Photon500_);
    obj.setPhoton600(HLT_Photon600_);
    obj.setMu300(HLT_Mu300_);
    obj.setMu350(HLT_Mu350_);
    obj.setMET250(HLT_MET250_);
    obj.setMET300(HLT_MET300_);
    obj.setMET600(HLT_MET600_);
    obj.setMET700(HLT_MET700_);
    obj.setPFMET300(HLT_PFMET300_);
    obj.setPFMET400(HLT_PFMET400_);
    obj.setPFMET500(HLT_PFMET500_);
    obj.setPFMET600(HLT_PFMET600_);
    obj.setEle250(HLT_Ele250_CaloIdVT_GsfTrkIdT_);
    obj.setEle300(HLT_Ele300_CaloIdVT_GsfTrkIdT_);
    obj.setHT2000(HLT_HT2000_);
    obj.setHT2500(HLT_HT2500_);
    obj.setIsoTrackHE(HLT_IsoTrackHE_);
    obj.setIsoTrackHB(HLT_IsoTrackHB_);
    obj.setHLTriggerFinalPath(HLTriggerFinalPath);
  
    return obj;
  }
  
  const vector<S>& SV(){
    if(SV_.size() > 0) return SV_;
    loadSs();
  	SV_.reserve(SV_dlen_->size());
    auto it_SV_dlen_ = SV_dlen_->cbegin();
    auto it_SV_dlenSig_ = SV_dlenSig_->cbegin();
    auto it_SV_pAngle_ = SV_pAngle_->cbegin();
    auto it_SV_chi2_ = SV_chi2_->cbegin();
    auto it_SV_eta_ = SV_eta_->cbegin();
    auto it_SV_mass_ = SV_mass_->cbegin();
    auto it_SV_ndof_ = SV_ndof_->cbegin();
    auto it_SV_phi_ = SV_phi_->cbegin();
    auto it_SV_pt_ = SV_pt_->cbegin();
    auto it_SV_x_ = SV_x_->cbegin();
    auto it_SV_y_ = SV_y_->cbegin();
    auto it_SV_z_ = SV_z_->cbegin();
    for(; it_SV_dlen_ != SV_dlen_->cend(); ){
      S obj;
      obj.setdlen(*it_SV_dlen_);
      obj.setdlenSig(*it_SV_dlenSig_);
      obj.setpAngle(*it_SV_pAngle_);
      obj.setchi2(*it_SV_chi2_);
      obj.setndof(*it_SV_ndof_);
      obj.setx(*it_SV_x_);
      obj.sety(*it_SV_y_);
      obj.setz(*it_SV_z_);
      obj.setLotentzVector(*it_SV_pt_, *it_SV_eta_, *it_SV_phi_, *it_SV_mass_);
      SV_.push_back( obj );
      ++it_SV_dlen_;
      ++it_SV_dlenSig_;
      ++it_SV_pAngle_;
      ++it_SV_chi2_;
      ++it_SV_eta_;
      ++it_SV_mass_;
      ++it_SV_ndof_;
      ++it_SV_phi_;
      ++it_SV_pt_;
      ++it_SV_x_;
      ++it_SV_y_;
      ++it_SV_z_;
    }
    return SV_;
  }
  
  const Met MET(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadMet();
  
    Met obj;
    obj.setMetUnclustEnUpDeltaX(MET_MetUnclustEnUpDeltaX_);
    obj.setMetUnclustEnUpDeltaY(MET_MetUnclustEnUpDeltaY_);
    obj.setcovXX(MET_covXX_);
    obj.setcovXY(MET_covXY_);
    obj.setcovYY(MET_covYY_);
    obj.setsignificance(MET_significance_);
    obj.setsumEt(MET_sumEt_);
    obj.setfiducialGenPhi(MET_fiducialGenPhi_);
    obj.setfiducialGenPt(MET_fiducialGenPt_);
  
    return obj;
  }
  
  const Genmet GenMET(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGenmet();
  
    Genmet obj;
    
  
    return obj;
  }
  
  const vector<Fatje>& FatJet(){
    if(FatJet_.size() > 0) return FatJet_;
    loadFatjes();
  	FatJet_.reserve(FatJet_area_->size());
    auto it_FatJet_area_ = FatJet_area_->cbegin();
    auto it_FatJet_btagCMVA_ = FatJet_btagCMVA_->cbegin();
    auto it_FatJet_btagCSVV2_ = FatJet_btagCSVV2_->cbegin();
    auto it_FatJet_btagDeepB_ = FatJet_btagDeepB_->cbegin();
    auto it_FatJet_btagHbb_ = FatJet_btagHbb_->cbegin();
    auto it_FatJet_eta_ = FatJet_eta_->cbegin();
    auto it_FatJet_mass_ = FatJet_mass_->cbegin();
    auto it_FatJet_msoftdrop_ = FatJet_msoftdrop_->cbegin();
    auto it_FatJet_msoftdrop_chs_ = FatJet_msoftdrop_chs_->cbegin();
    auto it_FatJet_n2b1_ = FatJet_n2b1_->cbegin();
    auto it_FatJet_n3b1_ = FatJet_n3b1_->cbegin();
    auto it_FatJet_phi_ = FatJet_phi_->cbegin();
    auto it_FatJet_pt_ = FatJet_pt_->cbegin();
    auto it_FatJet_tau1_ = FatJet_tau1_->cbegin();
    auto it_FatJet_tau2_ = FatJet_tau2_->cbegin();
    auto it_FatJet_tau3_ = FatJet_tau3_->cbegin();
    auto it_FatJet_tau4_ = FatJet_tau4_->cbegin();
    auto it_FatJet_jetId_ = FatJet_jetId_->cbegin();
    auto it_FatJet_subJetIdx1_ = FatJet_subJetIdx1_->cbegin();
    auto it_FatJet_subJetIdx2_ = FatJet_subJetIdx2_->cbegin();
    for(; it_FatJet_area_ != FatJet_area_->cend(); ){
      Fatje obj;
      obj.setarea(*it_FatJet_area_);
      obj.setbtagCMVA(*it_FatJet_btagCMVA_);
      obj.setbtagCSVV2(*it_FatJet_btagCSVV2_);
      obj.setbtagDeepB(*it_FatJet_btagDeepB_);
      obj.setbtagHbb(*it_FatJet_btagHbb_);
      obj.setmsoftdrop(*it_FatJet_msoftdrop_);
      obj.setmsoftdrop(*it_FatJet_msoftdrop_chs_);
      obj.setn2b1(*it_FatJet_n2b1_);
      obj.setn3b1(*it_FatJet_n3b1_);
      obj.settau1(*it_FatJet_tau1_);
      obj.settau2(*it_FatJet_tau2_);
      obj.settau3(*it_FatJet_tau3_);
      obj.settau4(*it_FatJet_tau4_);
      obj.setjetId(*it_FatJet_jetId_);
      obj.setsubJetIdx1(*it_FatJet_subJetIdx1_);
      obj.setsubJetIdx2(*it_FatJet_subJetIdx2_);
      obj.setLotentzVector(*it_FatJet_pt_, *it_FatJet_eta_, *it_FatJet_phi_, *it_FatJet_mass_);
      FatJet_.push_back( obj );
      ++it_FatJet_area_;
      ++it_FatJet_btagCMVA_;
      ++it_FatJet_btagCSVV2_;
      ++it_FatJet_btagDeepB_;
      ++it_FatJet_btagHbb_;
      ++it_FatJet_eta_;
      ++it_FatJet_mass_;
      ++it_FatJet_msoftdrop_;
      ++it_FatJet_msoftdrop_chs_;
      ++it_FatJet_n2b1_;
      ++it_FatJet_n3b1_;
      ++it_FatJet_phi_;
      ++it_FatJet_pt_;
      ++it_FatJet_tau1_;
      ++it_FatJet_tau2_;
      ++it_FatJet_tau3_;
      ++it_FatJet_tau4_;
      ++it_FatJet_jetId_;
      ++it_FatJet_subJetIdx1_;
      ++it_FatJet_subJetIdx2_;
    }
    return FatJet_;
  }
  
  const vector<Subje>& SubJet(){
    if(SubJet_.size() > 0) return SubJet_;
    loadSubjes();
  	SubJet_.reserve(SubJet_btagCMVA_->size());
    auto it_SubJet_btagCMVA_ = SubJet_btagCMVA_->cbegin();
    auto it_SubJet_btagCSVV2_ = SubJet_btagCSVV2_->cbegin();
    auto it_SubJet_btagDeepB_ = SubJet_btagDeepB_->cbegin();
    auto it_SubJet_eta_ = SubJet_eta_->cbegin();
    auto it_SubJet_mass_ = SubJet_mass_->cbegin();
    auto it_SubJet_n2b1_ = SubJet_n2b1_->cbegin();
    auto it_SubJet_n3b1_ = SubJet_n3b1_->cbegin();
    auto it_SubJet_phi_ = SubJet_phi_->cbegin();
    auto it_SubJet_pt_ = SubJet_pt_->cbegin();
    auto it_SubJet_tau1_ = SubJet_tau1_->cbegin();
    auto it_SubJet_tau2_ = SubJet_tau2_->cbegin();
    auto it_SubJet_tau3_ = SubJet_tau3_->cbegin();
    auto it_SubJet_tau4_ = SubJet_tau4_->cbegin();
    for(; it_SubJet_btagCMVA_ != SubJet_btagCMVA_->cend(); ){
      Subje obj;
      obj.setbtagCMVA(*it_SubJet_btagCMVA_);
      obj.setbtagCSVV2(*it_SubJet_btagCSVV2_);
      obj.setbtagDeepB(*it_SubJet_btagDeepB_);
      obj.setn2b1(*it_SubJet_n2b1_);
      obj.setn3b1(*it_SubJet_n3b1_);
      obj.settau1(*it_SubJet_tau1_);
      obj.settau2(*it_SubJet_tau2_);
      obj.settau3(*it_SubJet_tau3_);
      obj.settau4(*it_SubJet_tau4_);
      obj.setLotentzVector(*it_SubJet_pt_, *it_SubJet_eta_, *it_SubJet_phi_, *it_SubJet_mass_);
      SubJet_.push_back( obj );
      ++it_SubJet_btagCMVA_;
      ++it_SubJet_btagCSVV2_;
      ++it_SubJet_btagDeepB_;
      ++it_SubJet_eta_;
      ++it_SubJet_mass_;
      ++it_SubJet_n2b1_;
      ++it_SubJet_n3b1_;
      ++it_SubJet_phi_;
      ++it_SubJet_pt_;
      ++it_SubJet_tau1_;
      ++it_SubJet_tau2_;
      ++it_SubJet_tau3_;
      ++it_SubJet_tau4_;
    }
    return SubJet_;
  }
  
  const Flag Flag(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadFlag();
  
    Flag obj;
    obj.setHBHENoiseFilter(Flag_HBHENoiseFilter_);
    obj.setHBHENoiseIsoFilter(Flag_HBHENoiseIsoFilter_);
    obj.setCSCTightHaloFilter(Flag_CSCTightHaloFilter_);
    obj.setCSCTightHaloTrkMuUnvetoFilter(Flag_CSCTightHaloTrkMuUnvetoFilter_);
    obj.setCSCTightHalo2015Filter(Flag_CSCTightHalo2015Filter_);
    obj.setglobalTightHalo2016Filter(Flag_globalTightHalo2016Filter_);
    obj.setglobalSuperTightHalo2016Filter(Flag_globalSuperTightHalo2016Filter_);
    obj.setHcalStripHaloFilter(Flag_HcalStripHaloFilter_);
    obj.sethcalLaserEventFilter(Flag_hcalLaserEventFilter_);
    obj.setEcalDeadCellTriggerPrimitiveFilter(Flag_EcalDeadCellTriggerPrimitiveFilter_);
    obj.setEcalDeadCellBoundaryEnergyFilter(Flag_EcalDeadCellBoundaryEnergyFilter_);
    obj.setgoodVertices(Flag_goodVertices_);
    obj.seteeBadScFilter(Flag_eeBadScFilter_);
    obj.setecalLaserCorrFilter(Flag_ecalLaserCorrFilter_);
    obj.settrkPOGFilters(Flag_trkPOGFilters_);
    obj.setchargedHadronTrackResolutionFilter(Flag_chargedHadronTrackResolutionFilter_);
    obj.setmuonBadTrackFilter(Flag_muonBadTrackFilter_);
    obj.settrkPOG(Flag_trkPOG_manystripclus53X_);
    obj.settrkPOG(Flag_trkPOG_toomanystripclus53X_);
    obj.settrkPOG(Flag_trkPOG_logErrorTooManyClusters_);
    obj.setMETFilters(Flag_METFilters_);
  
    return obj;
  }
  
  const Pileup Pileup(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadPileup();
  
    Pileup obj;
    obj.setnTrueInt(Pileup_nTrueInt_);
    obj.setnPU(Pileup_nPU_);
    obj.setsumEOOT(Pileup_sumEOOT_);
    obj.setsumLOOT(Pileup_sumLOOT_);
  
    return obj;
  }
  
  const Lheweight LHEWeight(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadLheweight();
  
    Lheweight obj;
    obj.setoriginalXWGTUP(LHEWeight_originalXWGTUP_);
  
    return obj;
  }
  

private:
  TTree *tree_;
  Long64_t entries_;
  Long64_t current_entry_;
  Float_t CaloMET_phi_;
  Float_t CaloMET_pt_;
  Float_t CaloMET_sumEt_;
  Float_t *Electron_deltaEtaSC_;
  Float_t *Electron_dr03EcalRecHitSumEt_;
  Float_t *Electron_dr03HcalDepth1TowerSumEt_;
  Float_t *Electron_dr03TkSumPt_;
  Float_t *Electron_dxy_;
  Float_t *Electron_dxyErr_;
  Float_t *Electron_dz_;
  Float_t *Electron_dzErr_;
  Float_t *Electron_eCorr_;
  Float_t *Electron_eInvMinusPInv_;
  Float_t *Electron_energyErr_;
  Float_t *Electron_eta_;
  Float_t *Electron_hoe_;
  Float_t *Electron_ip3d_;
  Float_t *Electron_mass_;
  Float_t *Electron_miniPFRelIso_all_;
  Float_t *Electron_miniPFRelIso_chg_;
  Float_t *Electron_mvaSpring16GP_;
  Float_t *Electron_mvaSpring16HZZ_;
  Float_t *Electron_pfRelIso03_all_;
  Float_t *Electron_pfRelIso03_chg_;
  Float_t *Electron_phi_;
  Float_t *Electron_pt_;
  Float_t *Electron_r9_;
  Float_t *Electron_sieie_;
  Float_t *Electron_sip3d_;
  Float_t *Electron_mvaTTH_;
  Int_t *Electron_charge_;
  Int_t *Electron_cutBased_;
  Int_t *Electron_cutBased_HLTPreSel_;
  Int_t *Electron_jetIdx_;
  Int_t *Electron_pdgId_;
  Int_t *Electron_photonIdx_;
  Int_t *Electron_tightCharge_;
  Int_t *Electron_vidNestedWPBitmap_;
  Bool_t *Electron_convVeto_;
  Bool_t *Electron_cutBased_HEEP_;
  Bool_t *Electron_isPFcand_;
  UChar_t *Electron_lostHits_;
  Bool_t *Electron_mvaSpring16GP_WP80_;
  Bool_t *Electron_mvaSpring16GP_WP90_;
  Bool_t *Electron_mvaSpring16HZZ_WPL_;
  Float_t *FatJet_area_;
  Float_t *FatJet_btagCMVA_;
  Float_t *FatJet_btagCSVV2_;
  Float_t *FatJet_btagDeepB_;
  Float_t *FatJet_btagHbb_;
  Float_t *FatJet_eta_;
  Float_t *FatJet_mass_;
  Float_t *FatJet_msoftdrop_;
  Float_t *FatJet_msoftdrop_chs_;
  Float_t *FatJet_n2b1_;
  Float_t *FatJet_n3b1_;
  Float_t *FatJet_phi_;
  Float_t *FatJet_pt_;
  Float_t *FatJet_tau1_;
  Float_t *FatJet_tau2_;
  Float_t *FatJet_tau3_;
  Float_t *FatJet_tau4_;
  Int_t *FatJet_jetId_;
  Int_t *FatJet_subJetIdx1_;
  Int_t *FatJet_subJetIdx2_;
  Float_t *GenJetAK8_eta_;
  Float_t *GenJetAK8_mass_;
  Float_t *GenJetAK8_phi_;
  Float_t *GenJetAK8_pt_;
  Float_t *GenJet_eta_;
  Float_t *GenJet_mass_;
  Float_t *GenJet_phi_;
  Float_t *GenJet_pt_;
  Float_t *GenPart_eta_;
  Float_t *GenPart_mass_;
  Float_t *GenPart_phi_;
  Float_t *GenPart_pt_;
  Int_t *GenPart_genPartIdxMother_;
  Int_t *GenPart_pdgId_;
  Int_t *GenPart_status_;
  Int_t *GenPart_statusFlags_;
  Float_t Generator_binvar_;
  Float_t Generator_scalePDF_;
  Float_t Generator_weight_;
  Float_t Generator_x1_;
  Float_t Generator_x2_;
  Float_t Generator_xpdf1_;
  Float_t Generator_xpdf2_;
  Int_t Generator_id1_;
  Int_t Generator_id2_;
  Float_t *GenVisTau_eta_;
  Float_t *GenVisTau_mass_;
  Float_t *GenVisTau_phi_;
  Float_t *GenVisTau_pt_;
  Int_t *GenVisTau_charge_;
  Int_t *GenVisTau_genPartIdxMother_;
  Int_t *GenVisTau_status_;
  Float_t LHEWeight_originalXWGTUP_;
  Float_t *Jet_area_;
  Float_t *Jet_btagCMVA_;
  Float_t *Jet_btagCSVV2_;
  Float_t *Jet_btagDeepB_;
  Float_t *Jet_btagDeepC_;
  Float_t *Jet_chEmEF_;
  Float_t *Jet_chHEF_;
  Float_t *Jet_eta_;
  Float_t *Jet_mass_;
  Float_t *Jet_neEmEF_;
  Float_t *Jet_neHEF_;
  Float_t *Jet_phi_;
  Float_t *Jet_pt_;
  Float_t *Jet_qgl_;
  Float_t *Jet_rawFactor_;
  Float_t *Jet_bReg_;
  Int_t *Jet_electronIdx1_;
  Int_t *Jet_electronIdx2_;
  Int_t *Jet_jetId_;
  Int_t *Jet_muonIdx1_;
  Int_t *Jet_muonIdx2_;
  Int_t *Jet_nConstituents_;
  Int_t *Jet_nElectrons_;
  Int_t *Jet_nMuons_;
  Int_t *Jet_puId_;
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
  Float_t *Muon_dxy_;
  Float_t *Muon_dxyErr_;
  Float_t *Muon_dz_;
  Float_t *Muon_dzErr_;
  Float_t *Muon_eta_;
  Float_t *Muon_ip3d_;
  Float_t *Muon_mass_;
  Float_t *Muon_miniPFRelIso_all_;
  Float_t *Muon_miniPFRelIso_chg_;
  Float_t *Muon_pfRelIso03_all_;
  Float_t *Muon_pfRelIso03_chg_;
  Float_t *Muon_pfRelIso04_all_;
  Float_t *Muon_phi_;
  Float_t *Muon_pt_;
  Float_t *Muon_ptErr_;
  Float_t *Muon_segmentComp_;
  Float_t *Muon_sip3d_;
  Float_t *Muon_mvaTTH_;
  Int_t *Muon_charge_;
  Int_t *Muon_jetIdx_;
  Int_t *Muon_nStations_;
  Int_t *Muon_nTrackerLayers_;
  Int_t *Muon_pdgId_;
  Int_t *Muon_tightCharge_;
  UChar_t *Muon_highPtId_;
  Bool_t *Muon_isPFcand_;
  Bool_t *Muon_mediumId_;
  Bool_t *Muon_softId_;
  Bool_t *Muon_tightId_;
  Float_t *Photon_eCorr_;
  Float_t *Photon_energyErr_;
  Float_t *Photon_eta_;
  Float_t *Photon_hoe_;
  Float_t *Photon_mass_;
  Float_t *Photon_mvaID_;
  Float_t *Photon_pfRelIso03_all_;
  Float_t *Photon_pfRelIso03_chg_;
  Float_t *Photon_phi_;
  Float_t *Photon_pt_;
  Float_t *Photon_r9_;
  Float_t *Photon_sieie_;
  Int_t *Photon_charge_;
  Int_t *Photon_cutBased_;
  Int_t *Photon_electronIdx_;
  Int_t *Photon_jetIdx_;
  Int_t *Photon_pdgId_;
  Int_t *Photon_vidNestedWPBitmap_;
  Bool_t *Photon_electronVeto_;
  Bool_t *Photon_mvaID_WP80_;
  Bool_t *Photon_mvaID_WP90_;
  Bool_t *Photon_pixelSeed_;
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
  Float_t *GenDressedLepton_eta_;
  Float_t *GenDressedLepton_mass_;
  Float_t *GenDressedLepton_phi_;
  Float_t *GenDressedLepton_pt_;
  Int_t *GenDressedLepton_pdgId_;
  Float_t *SoftActivityJet_eta_;
  Float_t *SoftActivityJet_phi_;
  Float_t *SoftActivityJet_pt_;
  Float_t *SubJet_btagCMVA_;
  Float_t *SubJet_btagCSVV2_;
  Float_t *SubJet_btagDeepB_;
  Float_t *SubJet_eta_;
  Float_t *SubJet_mass_;
  Float_t *SubJet_n2b1_;
  Float_t *SubJet_n3b1_;
  Float_t *SubJet_phi_;
  Float_t *SubJet_pt_;
  Float_t *SubJet_tau1_;
  Float_t *SubJet_tau2_;
  Float_t *SubJet_tau3_;
  Float_t *SubJet_tau4_;
  Float_t *Tau_chargedIso_;
  Float_t *Tau_dxy_;
  Float_t *Tau_dz_;
  Float_t *Tau_eta_;
  Float_t *Tau_footprintCorr_;
  Float_t *Tau_leadTkDeltaEta_;
  Float_t *Tau_leadTkDeltaPhi_;
  Float_t *Tau_leadTkPtOverTauPt_;
  Float_t *Tau_mass_;
  Float_t *Tau_neutralIso_;
  Float_t *Tau_phi_;
  Float_t *Tau_photonsOutsideSignalCone_;
  Float_t *Tau_pt_;
  Float_t *Tau_puCorr_;
  Float_t *Tau_rawAntiEle_;
  Float_t *Tau_rawIso_;
  Float_t *Tau_rawMVAnewDM_;
  Float_t *Tau_rawMVAoldDM_;
  Float_t *Tau_rawMVAoldDMdR03_;
  Int_t *Tau_charge_;
  Int_t *Tau_decayMode_;
  Int_t *Tau_jetIdx_;
  Int_t *Tau_rawAntiEleCat_;
  UChar_t *Tau_idAntiEle_;
  UChar_t *Tau_idAntiMu_;
  Bool_t *Tau_idDecayMode_;
  Bool_t *Tau_idDecayModeNewDMs_;
  UChar_t *Tau_idMVAnewDM_;
  UChar_t *Tau_idMVAoldDM_;
  UChar_t *Tau_idMVAoldDMdR03_;
  Float_t TkMET_phi_;
  Float_t TkMET_pt_;
  Float_t TkMET_sumEt_;
  Float_t *TrigObj_pt_;
  Float_t *TrigObj_eta_;
  Float_t *TrigObj_phi_;
  Float_t *TrigObj_l1pt_;
  Float_t *TrigObj_l1pt_2_;
  Float_t *TrigObj_l2pt_;
  Int_t *TrigObj_id_;
  Int_t *TrigObj_l1iso_;
  Int_t *TrigObj_l1charge_;
  Int_t *TrigObj_filterBits_;
  Float_t *OtherPV_z_;
  Float_t PV_ndof_;
  Float_t PV_x_;
  Float_t PV_y_;
  Float_t PV_z_;
  Float_t PV_chi2_;
  Float_t PV_score_;
  Int_t PV_npvs_;
  Int_t PV_npvsGood_;
  Float_t *SV_dlen_;
  Float_t *SV_dlenSig_;
  Float_t *SV_pAngle_;
  Int_t *Electron_genPartIdx_;
  UChar_t *Electron_genPartFlav_;
  Int_t *GenJetAK8_partonFlavour_;
  UChar_t *GenJetAK8_hadronFlavour_;
  Int_t *GenJet_partonFlavour_;
  UChar_t *GenJet_hadronFlavour_;
  Int_t *Jet_genJetIdx_;
  Int_t *Jet_hadronFlavour_;
  Int_t *Jet_partonFlavour_;
  Int_t *Muon_genPartIdx_;
  UChar_t *Muon_genPartFlav_;
  Int_t *Photon_genPartIdx_;
  UChar_t *Photon_genPartFlav_;
  Float_t MET_fiducialGenPhi_;
  Float_t MET_fiducialGenPt_;
  UChar_t *Electron_cleanmask_;
  UChar_t *Jet_cleanmask_;
  UChar_t *Muon_cleanmask_;
  UChar_t *Photon_cleanmask_;
  UChar_t *Tau_cleanmask_;
  Float_t *SV_chi2_;
  Float_t *SV_eta_;
  Float_t *SV_mass_;
  Float_t *SV_ndof_;
  Float_t *SV_phi_;
  Float_t *SV_pt_;
  Float_t *SV_x_;
  Float_t *SV_y_;
  Float_t *SV_z_;
  Int_t *Tau_genPartIdx_;
  UChar_t *Tau_genPartFlav_;
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
  bool are_Jet_loaded_;
  vector<Je> Jet_;
  bool are_GenJetAK8_loaded_;
  vector<Genjetak> GenJetAK8_;
  bool are_GenVisTau_loaded_;
  vector<Genvista> GenVisTau_;
  bool are_CaloMET_loaded_;
  Calomet CaloMET_;
  bool are_GenDressedLepton_loaded_;
  vector<Gendressedlepto> GenDressedLepton_;
  bool are_PV_loaded_;
  Pv PV_;
  bool are_Generator_loaded_;
  Generator Generator_;
  bool are_TrigObj_loaded_;
  vector<Trigob> TrigObj_;
  bool are_Photon_loaded_;
  vector<Photo> Photon_;
  bool are_GenJet_loaded_;
  vector<Genje> GenJet_;
  bool are_RawMET_loaded_;
  Rawmet RawMET_;
  bool are_Electron_loaded_;
  vector<Electro> Electron_;
  bool are_SoftActivityJet_loaded_;
  Softactivityjet SoftActivityJet_;
  bool are_L1simulation_loaded_;
  L1Simulation L1simulation_;
  bool are_GenPart_loaded_;
  vector<Genpar> GenPart_;
  bool are_LHE_loaded_;
  Lhe LHE_;
  bool are_TkMET_loaded_;
  Tkmet TkMET_;
  bool are_Tau_loaded_;
  vector<Ta> Tau_;
  bool are_PuppiMET_loaded_;
  Puppimet PuppiMET_;
  bool are_Muon_loaded_;
  vector<Muo> Muon_;
  bool are_OtherPV_loaded_;
  vector<Otherp> OtherPV_;
  bool are_HLT_loaded_;
  Hlt HLT_;
  bool are_SV_loaded_;
  vector<S> SV_;
  bool are_MET_loaded_;
  Met MET_;
  bool are_GenMET_loaded_;
  Genmet GenMET_;
  bool are_FatJet_loaded_;
  vector<Fatje> FatJet_;
  bool are_SubJet_loaded_;
  vector<Subje> SubJet_;
  bool are_Flag_loaded_;
  Flag Flag_;
  bool are_Pileup_loaded_;
  Pileup Pileup_;
  bool are_LHEWeight_loaded_;
  Lheweight LHEWeight_;
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

