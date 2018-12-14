#ifndef URStreamer_h
#define URStreamer_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include <vector>
using namespace std;

class Jets: public TLorentzVector{
friend class URStreamer;
public:
//  Jets(const Float_t &i_area_,const Float_t &i_btagCMVA_,const Float_t &i_btagCSVV2_,const Float_t &i_btagDeepB_,const Float_t &i_btagDeepC_,const Float_t &i_btagDeepFlavB_,const Float_t &i_chEmEF_,const Float_t &i_chHEF_,const Float_t &i_muEF_,const Float_t &i_neEmEF_,const Float_t &i_neHEF_,const Float_t &i_qgl_,const Float_t &i_rawFactor_,const Float_t &i_bRegCorr_,const Float_t &i_bRegRes_,const Int_t &i_electronIdx1_,const Int_t &i_electronIdx2_,const Int_t &i_jetId_,const Int_t &i_muonIdx1_,const Int_t &i_muonIdx2_,const Int_t &i_nConstituents_,const Int_t &i_nElectrons_,const Int_t &i_nMuons_,const Int_t &i_puId_,const Int_t &i_genJetIdx_,const Int_t &i_hadronFlavour_,const Int_t &i_partonFlavour_,const UChar_t &i_cleanmask_):
//    
//  {}
  Jets():
    TLorentzVector(),
    area_(0),
    btagCMVA_(0),
    btagCSVV2_(0),
    btagDeepB_(0),
    btagDeepC_(0),
    btagDeepFlavB_(0),
    chEmEF_(0),
    chHEF_(0),
    muEF_(0),
    neEmEF_(0),
    neHEF_(0),
    qgl_(0),
    rawFactor_(0),
    bRegCorr_(0),
    bRegRes_(0),
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
  Float_t btagDeepFlavB() const {return btagDeepFlavB_;}
  Float_t chEmEF() const {return chEmEF_;}
  Float_t chHEF() const {return chHEF_;}
  Float_t muEF() const {return muEF_;}
  Float_t neEmEF() const {return neEmEF_;}
  Float_t neHEF() const {return neHEF_;}
  Float_t qgl() const {return qgl_;}
  Float_t rawFactor() const {return rawFactor_;}
  Float_t bRegCorr() const {return bRegCorr_;}
  Float_t bRegRes() const {return bRegRes_;}
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
  Float_t btagDeepFlavB_;
  Float_t chEmEF_;
  Float_t chHEF_;
  Float_t muEF_;
  Float_t neEmEF_;
  Float_t neHEF_;
  Float_t qgl_;
  Float_t rawFactor_;
  Float_t bRegCorr_;
  Float_t bRegRes_;
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
  void setbtagDeepFlavB(const Float_t value) {btagDeepFlavB_ = value;}
  void setchEmEF(const Float_t value) {chEmEF_ = value;}
  void setchHEF(const Float_t value) {chHEF_ = value;}
  void setmuEF(const Float_t value) {muEF_ = value;}
  void setneEmEF(const Float_t value) {neEmEF_ = value;}
  void setneHEF(const Float_t value) {neHEF_ = value;}
  void setqgl(const Float_t value) {qgl_ = value;}
  void setrawFactor(const Float_t value) {rawFactor_ = value;}
  void setbRegCorr(const Float_t value) {bRegCorr_ = value;}
  void setbRegRes(const Float_t value) {bRegRes_ = value;}
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
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Isotracks: public TLorentzVector{
friend class URStreamer;
public:
//  Isotracks(const Float_t &i_dxy_,const Float_t &i_dz_,const Float_t &i_pfRelIso03_all_,const Float_t &i_pfRelIso03_chg_,const Float_t &i_miniPFRelIso_all_,const Float_t &i_miniPFRelIso_chg_,const Int_t &i_pdgId_,const Bool_t &i_isHighPurityTrack_,const Bool_t &i_isPFcand_):
//    
//  {}
  Isotracks():
    TLorentzVector(),
    dxy_(0),
    dz_(0),
    pfRelIso03_all_(0),
    pfRelIso03_chg_(0),
    miniPFRelIso_all_(0),
    miniPFRelIso_chg_(0),
    pdgId_(0),
    isHighPurityTrack_(0),
    isPFcand_(0)
  {}
  Float_t dxy() const {return dxy_;}
  Float_t dz() const {return dz_;}
  Float_t pfRelIso03_all() const {return pfRelIso03_all_;}
  Float_t pfRelIso03_chg() const {return pfRelIso03_chg_;}
  Float_t miniPFRelIso_all() const {return miniPFRelIso_all_;}
  Float_t miniPFRelIso_chg() const {return miniPFRelIso_chg_;}
  Int_t pdgId() const {return pdgId_;}
  Bool_t isHighPurityTrack() const {return isHighPurityTrack_;}
  Bool_t isPFcand() const {return isPFcand_;}
private:
  Float_t dxy_;
  Float_t dz_;
  Float_t pfRelIso03_all_;
  Float_t pfRelIso03_chg_;
  Float_t miniPFRelIso_all_;
  Float_t miniPFRelIso_chg_;
  Int_t pdgId_;
  Bool_t isHighPurityTrack_;
  Bool_t isPFcand_;
  void setdxy(const Float_t value) {dxy_ = value;}
  void setdz(const Float_t value) {dz_ = value;}
  void setpfRelIso03_all(const Float_t value) {pfRelIso03_all_ = value;}
  void setpfRelIso03_chg(const Float_t value) {pfRelIso03_chg_ = value;}
  void setminiPFRelIso_all(const Float_t value) {miniPFRelIso_all_ = value;}
  void setminiPFRelIso_chg(const Float_t value) {miniPFRelIso_chg_ = value;}
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setisHighPurityTrack(const Bool_t value) {isHighPurityTrack_ = value;}
  void setisPFcand(const Bool_t value) {isPFcand_ = value;}
  void setLorentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Genjetak8s: public TLorentzVector{
friend class URStreamer;
public:
//  Genjetak8s(const Int_t &i_partonFlavour_,const UChar_t &i_hadronFlavour_):
//    
//  {}
  Genjetak8s():
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
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Genvistaus: public TLorentzVector{
friend class URStreamer;
public:
//  Genvistaus(const Int_t &i_charge_,const Int_t &i_genPartIdxMother_,const Int_t &i_status_):
//    
//  {}
  Genvistaus():
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
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
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
  void setLorentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Gendressedleptons: public TLorentzVector{
friend class URStreamer;
public:
//  Gendressedleptons(const Int_t &i_pdgId_):
//    
//  {}
  Gendressedleptons():
    TLorentzVector(),
    pdgId_(0)
  {}
  Int_t pdgId() const {return pdgId_;}
private:
  Int_t pdgId_;
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Lheparts: public TLorentzVector{
friend class URStreamer;
public:
//  Lheparts(const Int_t &i_pdgId_):
//    
//  {}
  Lheparts():
    TLorentzVector(),
    pdgId_(0)
  {}
  Int_t pdgId() const {return pdgId_;}
private:
  Int_t pdgId_;
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
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

class Trigobjs: public TLorentzVector{
friend class URStreamer;
public:
//  Trigobjs(const Float_t &i_l1pt_,const Float_t &i_l1pt_2_,const Float_t &i_l2pt_,const Int_t &i_id_,const Int_t &i_l1iso_,const Int_t &i_l1charge_,const Int_t &i_filterBits_):
//    
//  {}
  Trigobjs():
    TLorentzVector(),
    l1pt_(0),
    l1pt_2_(0),
    l2pt_(0),
    id_(0),
    l1iso_(0),
    l1charge_(0),
    filterBits_(0)
  {}
  Float_t l1pt() const {return l1pt_;}
  Float_t l1pt_2() const {return l1pt_2_;}
  Float_t l2pt() const {return l2pt_;}
  Int_t id() const {return id_;}
  Int_t l1iso() const {return l1iso_;}
  Int_t l1charge() const {return l1charge_;}
  Int_t filterBits() const {return filterBits_;}
private:
  Float_t l1pt_;
  Float_t l1pt_2_;
  Float_t l2pt_;
  Int_t id_;
  Int_t l1iso_;
  Int_t l1charge_;
  Int_t filterBits_;
  void setl1pt(const Float_t value) {l1pt_ = value;}
  void setl1pt_2(const Float_t value) {l1pt_2_ = value;}
  void setl2pt(const Float_t value) {l2pt_ = value;}
  void setid(const Int_t value) {id_ = value;}
  void setl1iso(const Int_t value) {l1iso_ = value;}
  void setl1charge(const Int_t value) {l1charge_ = value;}
  void setfilterBits(const Int_t value) {filterBits_ = value;}
  void setLorentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
};

class Photons: public TLorentzVector{
friend class URStreamer;
public:
//  Photons(const Float_t &i_energyErr_,const Float_t &i_hoe_,const Float_t &i_mvaID_,const Float_t &i_pfRelIso03_all_,const Float_t &i_pfRelIso03_chg_,const Float_t &i_r9_,const Float_t &i_sieie_,const Int_t &i_charge_,const Int_t &i_cutBasedBitmap_,const Int_t &i_electronIdx_,const Int_t &i_jetIdx_,const Int_t &i_pdgId_,const Int_t &i_vidNestedWPBitmap_,const Bool_t &i_electronVeto_,const Bool_t &i_isScEtaEB_,const Bool_t &i_isScEtaEE_,const Bool_t &i_mvaID_WP80_,const Bool_t &i_mvaID_WP90_,const Bool_t &i_pixelSeed_,const Int_t &i_genPartIdx_,const UChar_t &i_genPartFlav_,const UChar_t &i_cleanmask_):
//    
//  {}
  Photons():
    TLorentzVector(),
    energyErr_(0),
    hoe_(0),
    mvaID_(0),
    pfRelIso03_all_(0),
    pfRelIso03_chg_(0),
    r9_(0),
    sieie_(0),
    charge_(0),
    cutBasedBitmap_(0),
    electronIdx_(0),
    jetIdx_(0),
    pdgId_(0),
    vidNestedWPBitmap_(0),
    electronVeto_(0),
    isScEtaEB_(0),
    isScEtaEE_(0),
    mvaID_WP80_(0),
    mvaID_WP90_(0),
    pixelSeed_(0),
    genPartIdx_(0),
    genPartFlav_(0),
    cleanmask_(0)
  {}
  Float_t energyErr() const {return energyErr_;}
  Float_t hoe() const {return hoe_;}
  Float_t mvaID() const {return mvaID_;}
  Float_t pfRelIso03_all() const {return pfRelIso03_all_;}
  Float_t pfRelIso03_chg() const {return pfRelIso03_chg_;}
  Float_t r9() const {return r9_;}
  Float_t sieie() const {return sieie_;}
  Int_t charge() const {return charge_;}
  Int_t cutBasedBitmap() const {return cutBasedBitmap_;}
  Int_t electronIdx() const {return electronIdx_;}
  Int_t jetIdx() const {return jetIdx_;}
  Int_t pdgId() const {return pdgId_;}
  Int_t vidNestedWPBitmap() const {return vidNestedWPBitmap_;}
  Bool_t electronVeto() const {return electronVeto_;}
  Bool_t isScEtaEB() const {return isScEtaEB_;}
  Bool_t isScEtaEE() const {return isScEtaEE_;}
  Bool_t mvaID_WP80() const {return mvaID_WP80_;}
  Bool_t mvaID_WP90() const {return mvaID_WP90_;}
  Bool_t pixelSeed() const {return pixelSeed_;}
  Int_t genPartIdx() const {return genPartIdx_;}
  UChar_t genPartFlav() const {return genPartFlav_;}
  UChar_t cleanmask() const {return cleanmask_;}
private:
  Float_t energyErr_;
  Float_t hoe_;
  Float_t mvaID_;
  Float_t pfRelIso03_all_;
  Float_t pfRelIso03_chg_;
  Float_t r9_;
  Float_t sieie_;
  Int_t charge_;
  Int_t cutBasedBitmap_;
  Int_t electronIdx_;
  Int_t jetIdx_;
  Int_t pdgId_;
  Int_t vidNestedWPBitmap_;
  Bool_t electronVeto_;
  Bool_t isScEtaEB_;
  Bool_t isScEtaEE_;
  Bool_t mvaID_WP80_;
  Bool_t mvaID_WP90_;
  Bool_t pixelSeed_;
  Int_t genPartIdx_;
  UChar_t genPartFlav_;
  UChar_t cleanmask_;
  void setenergyErr(const Float_t value) {energyErr_ = value;}
  void sethoe(const Float_t value) {hoe_ = value;}
  void setmvaID(const Float_t value) {mvaID_ = value;}
  void setpfRelIso03_all(const Float_t value) {pfRelIso03_all_ = value;}
  void setpfRelIso03_chg(const Float_t value) {pfRelIso03_chg_ = value;}
  void setr9(const Float_t value) {r9_ = value;}
  void setsieie(const Float_t value) {sieie_ = value;}
  void setcharge(const Int_t value) {charge_ = value;}
  void setcutBasedBitmap(const Int_t value) {cutBasedBitmap_ = value;}
  void setelectronIdx(const Int_t value) {electronIdx_ = value;}
  void setjetIdx(const Int_t value) {jetIdx_ = value;}
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setvidNestedWPBitmap(const Int_t value) {vidNestedWPBitmap_ = value;}
  void setelectronVeto(const Bool_t value) {electronVeto_ = value;}
  void setisScEtaEB(const Bool_t value) {isScEtaEB_ = value;}
  void setisScEtaEE(const Bool_t value) {isScEtaEE_ = value;}
  void setmvaID_WP80(const Bool_t value) {mvaID_WP80_ = value;}
  void setmvaID_WP90(const Bool_t value) {mvaID_WP90_ = value;}
  void setpixelSeed(const Bool_t value) {pixelSeed_ = value;}
  void setgenPartIdx(const Int_t value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_t value) {genPartFlav_ = value;}
  void setcleanmask(const UChar_t value) {cleanmask_ = value;}
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Genjets: public TLorentzVector{
friend class URStreamer;
public:
//  Genjets(const Int_t &i_partonFlavour_,const UChar_t &i_hadronFlavour_,const Int_t &i_partonFlavour_,const UChar_t &i_hadronFlavour_):
//    
//  {}
  Genjets():
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
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
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
  void setLorentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Btagweight{
friend class URStreamer;
public:
//  Btagweight(const Float_t &i_CSVV2_,const Float_t &i_DeepCSVB_):
//    
//  {}
  Btagweight():
    CSVV2_(0),
    DeepCSVB_(0)
  {}
  Float_t CSVV2() const {return CSVV2_;}
  Float_t DeepCSVB() const {return DeepCSVB_;}
private:
  Float_t CSVV2_;
  Float_t DeepCSVB_;
  void setCSVV2(const Float_t value) {CSVV2_ = value;}
  void setDeepCSVB(const Float_t value) {DeepCSVB_ = value;}
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
  void setLorentzVector(float pt, float eta, float phi){SetPtEtaPhiM(pt, eta, phi, 0.);}
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

class Genparts: public TLorentzVector{
friend class URStreamer;
public:
//  Genparts(const Int_t &i_genPartIdxMother_,const Int_t &i_pdgId_,const Int_t &i_status_,const Int_t &i_statusFlags_):
//    
//  {}
  Genparts():
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
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Lhe: public TLorentzVector{
friend class URStreamer;
public:
//  Lhe(const Float_t &i_originalXWGTUP_,const vector<Float_t> &i_LHEPdfWeight_,const vector<Float_t> &i_LHEScaleWeight_,const Float_t &i_HT_,const Float_t &i_HTIncoming_,const Float_t &i_Vpt_,const UChar_t &i_Njets_,const UChar_t &i_Nb_,const UChar_t &i_Nc_,const UChar_t &i_Nuds_,const UChar_t &i_Nglu_,const UChar_t &i_NpNLO_,const UChar_t &i_NpLO_,const Int_t &i_pdgId_):
//    
//  {}
  Lhe():
    TLorentzVector(),
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
    NpLO_(0),
    pdgId_(0)
  {}
  vector<Float_t> LHEPdfWeight() const {return LHEPdfWeight_;}
  vector<Float_t> LHEScaleWeight() const {return LHEScaleWeight_;}
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
  Int_t pdgId() const {return pdgId_;}
private:
  vector<Float_t> LHEPdfWeight_;
  vector<Float_t> LHEScaleWeight_;
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
  Int_t pdgId_;
  void setLHEPdfWeight(const vector<Float_t> value) {LHEPdfWeight_ = value;}
  void setLHEScaleWeight(const vector<Float_t> value) {LHEScaleWeight_ = value;}
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
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
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
  void setLorentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Taus: public TLorentzVector{
friend class URStreamer;
public:
//  Taus(const Float_t &i_chargedIso_,const Float_t &i_dxy_,const Float_t &i_dz_,const Float_t &i_leadTkDeltaEta_,const Float_t &i_leadTkDeltaPhi_,const Float_t &i_leadTkPtOverTauPt_,const Float_t &i_neutralIso_,const Float_t &i_photonsOutsideSignalCone_,const Float_t &i_puCorr_,const Float_t &i_rawAntiEle_,const Float_t &i_rawIso_,const Float_t &i_rawIsodR03_,const Int_t &i_charge_,const Int_t &i_decayMode_,const Int_t &i_jetIdx_,const Int_t &i_rawAntiEleCat_,const UChar_t &i_idAntiEle_,const UChar_t &i_idAntiMu_,const Bool_t &i_idDecayMode_,const Bool_t &i_idDecayModeNewDMs_,const UChar_t &i_cleanmask_,const Int_t &i_genPartIdx_,const UChar_t &i_genPartFlav_):
//    
//  {}
  Taus():
    TLorentzVector(),
    chargedIso_(0),
    dxy_(0),
    dz_(0),
    leadTkDeltaEta_(0),
    leadTkDeltaPhi_(0),
    leadTkPtOverTauPt_(0),
    neutralIso_(0),
    photonsOutsideSignalCone_(0),
    puCorr_(0),
    rawAntiEle_(0),
    rawIso_(0),
    rawIsodR03_(0),
    charge_(0),
    decayMode_(0),
    jetIdx_(0),
    rawAntiEleCat_(0),
    idAntiEle_(0),
    idAntiMu_(0),
    idDecayMode_(0),
    idDecayModeNewDMs_(0),
    cleanmask_(0),
    genPartIdx_(0),
    genPartFlav_(0)
  {}
  Float_t chargedIso() const {return chargedIso_;}
  Float_t dxy() const {return dxy_;}
  Float_t dz() const {return dz_;}
  Float_t leadTkDeltaEta() const {return leadTkDeltaEta_;}
  Float_t leadTkDeltaPhi() const {return leadTkDeltaPhi_;}
  Float_t leadTkPtOverTauPt() const {return leadTkPtOverTauPt_;}
  Float_t neutralIso() const {return neutralIso_;}
  Float_t photonsOutsideSignalCone() const {return photonsOutsideSignalCone_;}
  Float_t puCorr() const {return puCorr_;}
  Float_t rawAntiEle() const {return rawAntiEle_;}
  Float_t rawIso() const {return rawIso_;}
  Float_t rawIsodR03() const {return rawIsodR03_;}
  Int_t charge() const {return charge_;}
  Int_t decayMode() const {return decayMode_;}
  Int_t jetIdx() const {return jetIdx_;}
  Int_t rawAntiEleCat() const {return rawAntiEleCat_;}
  UChar_t idAntiEle() const {return idAntiEle_;}
  UChar_t idAntiMu() const {return idAntiMu_;}
  Bool_t idDecayMode() const {return idDecayMode_;}
  Bool_t idDecayModeNewDMs() const {return idDecayModeNewDMs_;}
  UChar_t cleanmask() const {return cleanmask_;}
  Int_t genPartIdx() const {return genPartIdx_;}
  UChar_t genPartFlav() const {return genPartFlav_;}
private:
  Float_t chargedIso_;
  Float_t dxy_;
  Float_t dz_;
  Float_t leadTkDeltaEta_;
  Float_t leadTkDeltaPhi_;
  Float_t leadTkPtOverTauPt_;
  Float_t neutralIso_;
  Float_t photonsOutsideSignalCone_;
  Float_t puCorr_;
  Float_t rawAntiEle_;
  Float_t rawIso_;
  Float_t rawIsodR03_;
  Int_t charge_;
  Int_t decayMode_;
  Int_t jetIdx_;
  Int_t rawAntiEleCat_;
  UChar_t idAntiEle_;
  UChar_t idAntiMu_;
  Bool_t idDecayMode_;
  Bool_t idDecayModeNewDMs_;
  UChar_t cleanmask_;
  Int_t genPartIdx_;
  UChar_t genPartFlav_;
  void setchargedIso(const Float_t value) {chargedIso_ = value;}
  void setdxy(const Float_t value) {dxy_ = value;}
  void setdz(const Float_t value) {dz_ = value;}
  void setleadTkDeltaEta(const Float_t value) {leadTkDeltaEta_ = value;}
  void setleadTkDeltaPhi(const Float_t value) {leadTkDeltaPhi_ = value;}
  void setleadTkPtOverTauPt(const Float_t value) {leadTkPtOverTauPt_ = value;}
  void setneutralIso(const Float_t value) {neutralIso_ = value;}
  void setphotonsOutsideSignalCone(const Float_t value) {photonsOutsideSignalCone_ = value;}
  void setpuCorr(const Float_t value) {puCorr_ = value;}
  void setrawAntiEle(const Float_t value) {rawAntiEle_ = value;}
  void setrawIso(const Float_t value) {rawIso_ = value;}
  void setrawIsodR03(const Float_t value) {rawIsodR03_ = value;}
  void setcharge(const Int_t value) {charge_ = value;}
  void setdecayMode(const Int_t value) {decayMode_ = value;}
  void setjetIdx(const Int_t value) {jetIdx_ = value;}
  void setrawAntiEleCat(const Int_t value) {rawAntiEleCat_ = value;}
  void setidAntiEle(const UChar_t value) {idAntiEle_ = value;}
  void setidAntiMu(const UChar_t value) {idAntiMu_ = value;}
  void setidDecayMode(const Bool_t value) {idDecayMode_ = value;}
  void setidDecayModeNewDMs(const Bool_t value) {idDecayModeNewDMs_ = value;}
  void setcleanmask(const UChar_t value) {cleanmask_ = value;}
  void setgenPartIdx(const Int_t value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_t value) {genPartFlav_ = value;}
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Flag{
friend class URStreamer;
public:
//  Flag(const Bool_t &i_HBHENoiseFilter_,const Bool_t &i_HBHENoiseIsoFilter_,const Bool_t &i_CSCTightHaloFilter_,const Bool_t &i_CSCTightHaloTrkMuUnvetoFilter_,const Bool_t &i_CSCTightHalo2015Filter_,const Bool_t &i_globalTightHalo2016Filter_,const Bool_t &i_globalSuperTightHalo2016Filter_,const Bool_t &i_HcalStripHaloFilter_,const Bool_t &i_hcalLaserEventFilter_,const Bool_t &i_EcalDeadCellTriggerPrimitiveFilter_,const Bool_t &i_EcalDeadCellBoundaryEnergyFilter_,const Bool_t &i_ecalBadCalibFilter_,const Bool_t &i_goodVertices_,const Bool_t &i_eeBadScFilter_,const Bool_t &i_ecalLaserCorrFilter_,const Bool_t &i_trkPOGFilters_,const Bool_t &i_chargedHadronTrackResolutionFilter_,const Bool_t &i_muonBadTrackFilter_,const Bool_t &i_BadChargedCandidateFilter_,const Bool_t &i_BadPFMuonFilter_,const Bool_t &i_BadChargedCandidateSummer16Filter_,const Bool_t &i_BadPFMuonSummer16Filter_,const Bool_t &i_trkPOG_manystripclus53X_,const Bool_t &i_trkPOG_toomanystripclus53X_,const Bool_t &i_trkPOG_logErrorTooManyClusters_,const Bool_t &i_METFilters_):
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
    ecalBadCalibFilter_(0),
    goodVertices_(0),
    eeBadScFilter_(0),
    ecalLaserCorrFilter_(0),
    trkPOGFilters_(0),
    chargedHadronTrackResolutionFilter_(0),
    muonBadTrackFilter_(0),
    BadChargedCandidateFilter_(0),
    BadPFMuonFilter_(0),
    BadChargedCandidateSummer16Filter_(0),
    BadPFMuonSummer16Filter_(0),
    trkPOG_manystripclus53X_(0),
    trkPOG_toomanystripclus53X_(0),
    trkPOG_logErrorTooManyClusters_(0),
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
  Bool_t ecalBadCalibFilter() const {return ecalBadCalibFilter_;}
  Bool_t goodVertices() const {return goodVertices_;}
  Bool_t eeBadScFilter() const {return eeBadScFilter_;}
  Bool_t ecalLaserCorrFilter() const {return ecalLaserCorrFilter_;}
  Bool_t trkPOGFilters() const {return trkPOGFilters_;}
  Bool_t chargedHadronTrackResolutionFilter() const {return chargedHadronTrackResolutionFilter_;}
  Bool_t muonBadTrackFilter() const {return muonBadTrackFilter_;}
  Bool_t BadChargedCandidateFilter() const {return BadChargedCandidateFilter_;}
  Bool_t BadPFMuonFilter() const {return BadPFMuonFilter_;}
  Bool_t BadChargedCandidateSummer16Filter() const {return BadChargedCandidateSummer16Filter_;}
  Bool_t BadPFMuonSummer16Filter() const {return BadPFMuonSummer16Filter_;}
  Bool_t trkPOG_manystripclus53X() const {return trkPOG_manystripclus53X_;}
  Bool_t trkPOG_toomanystripclus53X() const {return trkPOG_toomanystripclus53X_;}
  Bool_t trkPOG_logErrorTooManyClusters() const {return trkPOG_logErrorTooManyClusters_;}
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
  Bool_t ecalBadCalibFilter_;
  Bool_t goodVertices_;
  Bool_t eeBadScFilter_;
  Bool_t ecalLaserCorrFilter_;
  Bool_t trkPOGFilters_;
  Bool_t chargedHadronTrackResolutionFilter_;
  Bool_t muonBadTrackFilter_;
  Bool_t BadChargedCandidateFilter_;
  Bool_t BadPFMuonFilter_;
  Bool_t BadChargedCandidateSummer16Filter_;
  Bool_t BadPFMuonSummer16Filter_;
  Bool_t trkPOG_manystripclus53X_;
  Bool_t trkPOG_toomanystripclus53X_;
  Bool_t trkPOG_logErrorTooManyClusters_;
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
  void setecalBadCalibFilter(const Bool_t value) {ecalBadCalibFilter_ = value;}
  void setgoodVertices(const Bool_t value) {goodVertices_ = value;}
  void seteeBadScFilter(const Bool_t value) {eeBadScFilter_ = value;}
  void setecalLaserCorrFilter(const Bool_t value) {ecalLaserCorrFilter_ = value;}
  void settrkPOGFilters(const Bool_t value) {trkPOGFilters_ = value;}
  void setchargedHadronTrackResolutionFilter(const Bool_t value) {chargedHadronTrackResolutionFilter_ = value;}
  void setmuonBadTrackFilter(const Bool_t value) {muonBadTrackFilter_ = value;}
  void setBadChargedCandidateFilter(const Bool_t value) {BadChargedCandidateFilter_ = value;}
  void setBadPFMuonFilter(const Bool_t value) {BadPFMuonFilter_ = value;}
  void setBadChargedCandidateSummer16Filter(const Bool_t value) {BadChargedCandidateSummer16Filter_ = value;}
  void setBadPFMuonSummer16Filter(const Bool_t value) {BadPFMuonSummer16Filter_ = value;}
  void settrkPOG_manystripclus53X(const Bool_t value) {trkPOG_manystripclus53X_ = value;}
  void settrkPOG_toomanystripclus53X(const Bool_t value) {trkPOG_toomanystripclus53X_ = value;}
  void settrkPOG_logErrorTooManyClusters(const Bool_t value) {trkPOG_logErrorTooManyClusters_ = value;}
  void setMETFilters(const Bool_t value) {METFilters_ = value;}
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
  void setLorentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Muons: public TLorentzVector{
friend class URStreamer;
public:
//  Muons(const Float_t &i_dxy_,const Float_t &i_dxyErr_,const Float_t &i_dz_,const Float_t &i_dzErr_,const Float_t &i_ip3d_,const Float_t &i_miniPFRelIso_all_,const Float_t &i_miniPFRelIso_chg_,const Float_t &i_pfRelIso03_all_,const Float_t &i_pfRelIso03_chg_,const Float_t &i_pfRelIso04_all_,const Float_t &i_ptErr_,const Float_t &i_segmentComp_,const Float_t &i_sip3d_,const Float_t &i_mvaTTH_,const Int_t &i_charge_,const Int_t &i_jetIdx_,const Int_t &i_nStations_,const Int_t &i_nTrackerLayers_,const Int_t &i_pdgId_,const Int_t &i_tightCharge_,const UChar_t &i_highPtId_,const Bool_t &i_isPFcand_,const Bool_t &i_mediumId_,const Bool_t &i_softId_,const Bool_t &i_tightId_,const Int_t &i_genPartIdx_,const UChar_t &i_genPartFlav_,const UChar_t &i_cleanmask_):
//    
//  {}
  Muons():
    TLorentzVector(),
    dxy_(0),
    dxyErr_(0),
    dz_(0),
    dzErr_(0),
    ip3d_(0),
    miniPFRelIso_all_(0),
    miniPFRelIso_chg_(0),
    pfRelIso03_all_(0),
    pfRelIso03_chg_(0),
    pfRelIso04_all_(0),
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
  Float_t miniPFRelIso_all() const {return miniPFRelIso_all_;}
  Float_t miniPFRelIso_chg() const {return miniPFRelIso_chg_;}
  Float_t pfRelIso03_all() const {return pfRelIso03_all_;}
  Float_t pfRelIso03_chg() const {return pfRelIso03_chg_;}
  Float_t pfRelIso04_all() const {return pfRelIso04_all_;}
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
  Float_t miniPFRelIso_all_;
  Float_t miniPFRelIso_chg_;
  Float_t pfRelIso03_all_;
  Float_t pfRelIso03_chg_;
  Float_t pfRelIso04_all_;
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
  void setminiPFRelIso_all(const Float_t value) {miniPFRelIso_all_ = value;}
  void setminiPFRelIso_chg(const Float_t value) {miniPFRelIso_chg_ = value;}
  void setpfRelIso03_all(const Float_t value) {pfRelIso03_all_ = value;}
  void setpfRelIso03_chg(const Float_t value) {pfRelIso03_chg_ = value;}
  void setpfRelIso04_all(const Float_t value) {pfRelIso04_all_ = value;}
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
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class L1Reco{
friend class URStreamer;
public:
//  L1Reco(const Bool_t &i_step_):
//    
//  {}
  L1Reco():
    step_(0)
  {}
  Bool_t step() const {return step_;}
private:
  Bool_t step_;
  void setstep(const Bool_t value) {step_ = value;}
};

class Otherpvs{
friend class URStreamer;
public:
//  Otherpvs(const Float_t &i_z_):
//    
//  {}
  Otherpvs():
    z_(0)
  {}
  Float_t z() const {return z_;}
private:
  Float_t z_;
  void setz(const Float_t value) {z_ = value;}
};

class Hlt{
friend class URStreamer;
public:
//  Hlt(const Bool_t &i_HLTriggerFirstPath_,const Bool_t &i_AK8PFJet360_TrimMass30_,const Bool_t &i_AK8PFJet380_TrimMass30_,const Bool_t &i_AK8PFJet400_TrimMass30_,const Bool_t &i_AK8PFJet420_TrimMass30_,const Bool_t &i_AK8PFHT750_TrimMass50_,const Bool_t &i_AK8PFHT800_TrimMass50_,const Bool_t &i_AK8PFHT850_TrimMass50_,const Bool_t &i_AK8PFHT900_TrimMass50_,const Bool_t &i_CaloJet500_NoJetID_,const Bool_t &i_CaloJet550_NoJetID_,const Bool_t &i_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_,const Bool_t &i_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_,const Bool_t &i_Trimuon5_3p5_2_Upsilon_Muon_,const Bool_t &i_TrimuonOpen_5_3p5_2_Upsilon_Muon_,const Bool_t &i_DoubleEle25_CaloIdL_MW_,const Bool_t &i_DoubleEle27_CaloIdL_MW_,const Bool_t &i_DoubleEle33_CaloIdL_MW_,const Bool_t &i_DoubleEle24_eta2p1_WPTight_Gsf_,const Bool_t &i_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_,const Bool_t &i_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_,const Bool_t &i_Ele27_Ele37_CaloIdL_MW_,const Bool_t &i_Mu27_Ele37_CaloIdL_MW_,const Bool_t &i_Mu37_Ele27_CaloIdL_MW_,const Bool_t &i_Mu37_TkMu27_,const Bool_t &i_DoubleMu4_3_Bs_,const Bool_t &i_DoubleMu4_3_Jpsi_,const Bool_t &i_DoubleMu4_JpsiTrk_Displaced_,const Bool_t &i_DoubleMu4_LowMassNonResonantTrk_Displaced_,const Bool_t &i_DoubleMu3_Trk_Tau3mu_,const Bool_t &i_DoubleMu3_TkMu_DsTau3Mu_,const Bool_t &i_DoubleMu4_PsiPrimeTrk_Displaced_,const Bool_t &i_DoubleMu4_Mass3p8_DZ_PFHT350_,const Bool_t &i_Mu3_PFJet40_,const Bool_t &i_Mu7p5_L2Mu2_Jpsi_,const Bool_t &i_Mu7p5_L2Mu2_Upsilon_,const Bool_t &i_Mu7p5_Track2_Jpsi_,const Bool_t &i_Mu7p5_Track3p5_Jpsi_,const Bool_t &i_Mu7p5_Track7_Jpsi_,const Bool_t &i_Mu7p5_Track2_Upsilon_,const Bool_t &i_Mu7p5_Track3p5_Upsilon_,const Bool_t &i_Mu7p5_Track7_Upsilon_,const Bool_t &i_Mu3_L1SingleMu5orSingleMu7_,const Bool_t &i_DoublePhoton33_CaloIdL_,const Bool_t &i_DoublePhoton70_,const Bool_t &i_DoublePhoton85_,const Bool_t &i_Ele20_WPTight_Gsf_,const Bool_t &i_Ele15_WPLoose_Gsf_,const Bool_t &i_Ele17_WPLoose_Gsf_,const Bool_t &i_Ele20_WPLoose_Gsf_,const Bool_t &i_Ele20_eta2p1_WPLoose_Gsf_,const Bool_t &i_DiEle27_WPTightCaloOnly_L1DoubleEG_,const Bool_t &i_Ele27_WPTight_Gsf_,const Bool_t &i_Ele28_WPTight_Gsf_,const Bool_t &i_Ele30_WPTight_Gsf_,const Bool_t &i_Ele32_WPTight_Gsf_,const Bool_t &i_Ele35_WPTight_Gsf_,const Bool_t &i_Ele35_WPTight_Gsf_L1EGMT_,const Bool_t &i_Ele38_WPTight_Gsf_,const Bool_t &i_Ele40_WPTight_Gsf_,const Bool_t &i_Ele32_WPTight_Gsf_L1DoubleEG_,const Bool_t &i_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_,const Bool_t &i_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_,const Bool_t &i_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_,const Bool_t &i_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_,const Bool_t &i_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_,const Bool_t &i_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_,const Bool_t &i_HT450_Beamspot_,const Bool_t &i_HT300_Beamspot_,const Bool_t &i_ZeroBias_Beamspot_,const Bool_t &i_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_,const Bool_t &i_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_,const Bool_t &i_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_,const Bool_t &i_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_,const Bool_t &i_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_,const Bool_t &i_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_,const Bool_t &i_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_,const Bool_t &i_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_,const Bool_t &i_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_,const Bool_t &i_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_,const Bool_t &i_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_,const Bool_t &i_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_,const Bool_t &i_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_,const Bool_t &i_IsoMu20_,const Bool_t &i_IsoMu24_,const Bool_t &i_IsoMu24_eta2p1_,const Bool_t &i_IsoMu27_,const Bool_t &i_IsoMu30_,const Bool_t &i_UncorrectedJetE30_NoBPTX_,const Bool_t &i_UncorrectedJetE30_NoBPTX3BX_,const Bool_t &i_UncorrectedJetE60_NoBPTX3BX_,const Bool_t &i_UncorrectedJetE70_NoBPTX3BX_,const Bool_t &i_L1SingleMu18_,const Bool_t &i_L1SingleMu25_,const Bool_t &i_L2Mu10_,const Bool_t &i_L2Mu10_NoVertex_NoBPTX3BX_,const Bool_t &i_L2Mu10_NoVertex_NoBPTX_,const Bool_t &i_L2Mu45_NoVertex_3Sta_NoBPTX3BX_,const Bool_t &i_L2Mu40_NoVertex_3Sta_NoBPTX3BX_,const Bool_t &i_L2Mu50_,const Bool_t &i_L2Mu23NoVtx_2Cha_,const Bool_t &i_L2Mu23NoVtx_2Cha_CosmicSeed_,const Bool_t &i_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_,const Bool_t &i_DoubleL2Mu30NoVtx_2Cha_Eta2p4_,const Bool_t &i_DoubleL2Mu50_,const Bool_t &i_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_,const Bool_t &i_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_,const Bool_t &i_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_,const Bool_t &i_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_,const Bool_t &i_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_,const Bool_t &i_DoubleL2Mu23NoVtx_2Cha_,const Bool_t &i_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_,const Bool_t &i_DoubleL2Mu25NoVtx_2Cha_,const Bool_t &i_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_,const Bool_t &i_DoubleL2Mu25NoVtx_2Cha_Eta2p4_,const Bool_t &i_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_,const Bool_t &i_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_,const Bool_t &i_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_,const Bool_t &i_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_,const Bool_t &i_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_,const Bool_t &i_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_,const Bool_t &i_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_,const Bool_t &i_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_,const Bool_t &i_Mu25_TkMu0_Onia_,const Bool_t &i_Mu30_TkMu0_Psi_,const Bool_t &i_Mu30_TkMu0_Upsilon_,const Bool_t &i_Mu20_TkMu0_Phi_,const Bool_t &i_Mu25_TkMu0_Phi_,const Bool_t &i_Mu12_,const Bool_t &i_Mu15_,const Bool_t &i_Mu20_,const Bool_t &i_Mu27_,const Bool_t &i_Mu50_,const Bool_t &i_Mu55_,const Bool_t &i_OldMu100_,const Bool_t &i_TkMu100_,const Bool_t &i_DiPFJetAve40_,const Bool_t &i_DiPFJetAve60_,const Bool_t &i_DiPFJetAve80_,const Bool_t &i_DiPFJetAve140_,const Bool_t &i_DiPFJetAve200_,const Bool_t &i_DiPFJetAve260_,const Bool_t &i_DiPFJetAve320_,const Bool_t &i_DiPFJetAve400_,const Bool_t &i_DiPFJetAve500_,const Bool_t &i_DiPFJetAve60_HFJEC_,const Bool_t &i_DiPFJetAve80_HFJEC_,const Bool_t &i_DiPFJetAve100_HFJEC_,const Bool_t &i_DiPFJetAve160_HFJEC_,const Bool_t &i_DiPFJetAve220_HFJEC_,const Bool_t &i_DiPFJetAve300_HFJEC_,const Bool_t &i_AK8PFJet15_,const Bool_t &i_AK8PFJet25_,const Bool_t &i_AK8PFJet40_,const Bool_t &i_AK8PFJet60_,const Bool_t &i_AK8PFJet80_,const Bool_t &i_AK8PFJet140_,const Bool_t &i_AK8PFJet200_,const Bool_t &i_AK8PFJet260_,const Bool_t &i_AK8PFJet320_,const Bool_t &i_AK8PFJet400_,const Bool_t &i_AK8PFJet450_,const Bool_t &i_AK8PFJet500_,const Bool_t &i_AK8PFJet550_,const Bool_t &i_PFJet15_,const Bool_t &i_PFJet25_,const Bool_t &i_PFJet40_,const Bool_t &i_PFJet60_,const Bool_t &i_PFJet80_,const Bool_t &i_PFJet140_,const Bool_t &i_PFJet200_,const Bool_t &i_PFJet260_,const Bool_t &i_PFJet320_,const Bool_t &i_PFJet400_,const Bool_t &i_PFJet450_,const Bool_t &i_PFJet500_,const Bool_t &i_PFJet550_,const Bool_t &i_PFJetFwd15_,const Bool_t &i_PFJetFwd25_,const Bool_t &i_PFJetFwd40_,const Bool_t &i_PFJetFwd60_,const Bool_t &i_PFJetFwd80_,const Bool_t &i_PFJetFwd140_,const Bool_t &i_PFJetFwd200_,const Bool_t &i_PFJetFwd260_,const Bool_t &i_PFJetFwd320_,const Bool_t &i_PFJetFwd400_,const Bool_t &i_PFJetFwd450_,const Bool_t &i_PFJetFwd500_,const Bool_t &i_AK8PFJetFwd15_,const Bool_t &i_AK8PFJetFwd25_,const Bool_t &i_AK8PFJetFwd40_,const Bool_t &i_AK8PFJetFwd60_,const Bool_t &i_AK8PFJetFwd80_,const Bool_t &i_AK8PFJetFwd140_,const Bool_t &i_AK8PFJetFwd200_,const Bool_t &i_AK8PFJetFwd260_,const Bool_t &i_AK8PFJetFwd320_,const Bool_t &i_AK8PFJetFwd400_,const Bool_t &i_AK8PFJetFwd450_,const Bool_t &i_AK8PFJetFwd500_,const Bool_t &i_PFHT180_,const Bool_t &i_PFHT250_,const Bool_t &i_PFHT370_,const Bool_t &i_PFHT430_,const Bool_t &i_PFHT510_,const Bool_t &i_PFHT590_,const Bool_t &i_PFHT680_,const Bool_t &i_PFHT780_,const Bool_t &i_PFHT890_,const Bool_t &i_PFHT1050_,const Bool_t &i_PFHT500_PFMET100_PFMHT100_IDTight_,const Bool_t &i_PFHT500_PFMET110_PFMHT110_IDTight_,const Bool_t &i_PFHT700_PFMET85_PFMHT85_IDTight_,const Bool_t &i_PFHT700_PFMET95_PFMHT95_IDTight_,const Bool_t &i_PFHT800_PFMET75_PFMHT75_IDTight_,const Bool_t &i_PFHT800_PFMET85_PFMHT85_IDTight_,const Bool_t &i_PFMET110_PFMHT110_IDTight_,const Bool_t &i_PFMET120_PFMHT120_IDTight_,const Bool_t &i_PFMET130_PFMHT130_IDTight_,const Bool_t &i_PFMET140_PFMHT140_IDTight_,const Bool_t &i_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_,const Bool_t &i_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_,const Bool_t &i_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_,const Bool_t &i_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_,const Bool_t &i_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_,const Bool_t &i_PFMET120_PFMHT120_IDTight_PFHT60_,const Bool_t &i_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_,const Bool_t &i_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_,const Bool_t &i_PFMETTypeOne110_PFMHT110_IDTight_,const Bool_t &i_PFMETTypeOne120_PFMHT120_IDTight_,const Bool_t &i_PFMETTypeOne130_PFMHT130_IDTight_,const Bool_t &i_PFMETTypeOne140_PFMHT140_IDTight_,const Bool_t &i_PFMETNoMu110_PFMHTNoMu110_IDTight_,const Bool_t &i_PFMETNoMu120_PFMHTNoMu120_IDTight_,const Bool_t &i_PFMETNoMu130_PFMHTNoMu130_IDTight_,const Bool_t &i_PFMETNoMu140_PFMHTNoMu140_IDTight_,const Bool_t &i_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_,const Bool_t &i_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_,const Bool_t &i_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_,const Bool_t &i_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_,const Bool_t &i_L1ETMHadSeeds_,const Bool_t &i_CaloMHT90_,const Bool_t &i_CaloMET80_NotCleaned_,const Bool_t &i_CaloMET90_NotCleaned_,const Bool_t &i_CaloMET100_NotCleaned_,const Bool_t &i_CaloMET110_NotCleaned_,const Bool_t &i_CaloMET250_NotCleaned_,const Bool_t &i_CaloMET70_HBHECleaned_,const Bool_t &i_CaloMET80_HBHECleaned_,const Bool_t &i_CaloMET90_HBHECleaned_,const Bool_t &i_CaloMET100_HBHECleaned_,const Bool_t &i_CaloMET250_HBHECleaned_,const Bool_t &i_CaloMET300_HBHECleaned_,const Bool_t &i_CaloMET350_HBHECleaned_,const Bool_t &i_PFMET200_NotCleaned_,const Bool_t &i_PFMET200_HBHECleaned_,const Bool_t &i_PFMET250_HBHECleaned_,const Bool_t &i_PFMET300_HBHECleaned_,const Bool_t &i_PFMET200_HBHE_BeamHaloCleaned_,const Bool_t &i_PFMETTypeOne200_HBHE_BeamHaloCleaned_,const Bool_t &i_MET105_IsoTrk50_,const Bool_t &i_MET120_IsoTrk50_,const Bool_t &i_SingleJet30_Mu12_SinglePFJet40_,const Bool_t &i_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_,const Bool_t &i_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_,const Bool_t &i_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_,const Bool_t &i_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_,const Bool_t &i_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_,const Bool_t &i_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_,const Bool_t &i_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_,const Bool_t &i_DoublePFJets40_CaloBTagDeepCSV_p71_,const Bool_t &i_DoublePFJets100_CaloBTagDeepCSV_p71_,const Bool_t &i_DoublePFJets200_CaloBTagDeepCSV_p71_,const Bool_t &i_DoublePFJets350_CaloBTagDeepCSV_p71_,const Bool_t &i_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_,const Bool_t &i_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_,const Bool_t &i_Photon300_NoHE_,const Bool_t &i_Mu8_TrkIsoVVL_,const Bool_t &i_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_,const Bool_t &i_Mu8_DiEle12_CaloIdL_TrackIdL_,const Bool_t &i_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_,const Bool_t &i_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_,const Bool_t &i_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_,const Bool_t &i_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_,const Bool_t &i_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_,const Bool_t &i_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_,const Bool_t &i_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_,const Bool_t &i_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_,const Bool_t &i_Mu17_TrkIsoVVL_,const Bool_t &i_Mu19_TrkIsoVVL_,const Bool_t &i_BTagMu_AK4DiJet20_Mu5_,const Bool_t &i_BTagMu_AK4DiJet40_Mu5_,const Bool_t &i_BTagMu_AK4DiJet70_Mu5_,const Bool_t &i_BTagMu_AK4DiJet110_Mu5_,const Bool_t &i_BTagMu_AK4DiJet170_Mu5_,const Bool_t &i_BTagMu_AK4Jet300_Mu5_,const Bool_t &i_BTagMu_AK8DiJet170_Mu5_,const Bool_t &i_BTagMu_AK8Jet170_DoubleMu5_,const Bool_t &i_BTagMu_AK8Jet300_Mu5_,const Bool_t &i_BTagMu_AK4DiJet20_Mu5_noalgo_,const Bool_t &i_BTagMu_AK4DiJet40_Mu5_noalgo_,const Bool_t &i_BTagMu_AK4DiJet70_Mu5_noalgo_,const Bool_t &i_BTagMu_AK4DiJet110_Mu5_noalgo_,const Bool_t &i_BTagMu_AK4DiJet170_Mu5_noalgo_,const Bool_t &i_BTagMu_AK4Jet300_Mu5_noalgo_,const Bool_t &i_BTagMu_AK8DiJet170_Mu5_noalgo_,const Bool_t &i_BTagMu_AK8Jet170_DoubleMu5_noalgo_,const Bool_t &i_BTagMu_AK8Jet300_Mu5_noalgo_,const Bool_t &i_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_,const Bool_t &i_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_,const Bool_t &i_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_,const Bool_t &i_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_,const Bool_t &i_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_,const Bool_t &i_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_,const Bool_t &i_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_,const Bool_t &i_Mu12_DoublePhoton20_,const Bool_t &i_TriplePhoton_20_20_20_CaloIdLV2_,const Bool_t &i_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_,const Bool_t &i_TriplePhoton_30_30_10_CaloIdLV2_,const Bool_t &i_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_,const Bool_t &i_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_,const Bool_t &i_Photon20_,const Bool_t &i_Photon33_,const Bool_t &i_Photon50_,const Bool_t &i_Photon75_,const Bool_t &i_Photon90_,const Bool_t &i_Photon120_,const Bool_t &i_Photon150_,const Bool_t &i_Photon175_,const Bool_t &i_Photon200_,const Bool_t &i_Photon100EB_TightID_TightIso_,const Bool_t &i_Photon110EB_TightID_TightIso_,const Bool_t &i_Photon120EB_TightID_TightIso_,const Bool_t &i_Photon100EBHE10_,const Bool_t &i_Photon100EEHE10_,const Bool_t &i_Photon100EE_TightID_TightIso_,const Bool_t &i_Photon50_R9Id90_HE10_IsoM_,const Bool_t &i_Photon75_R9Id90_HE10_IsoM_,const Bool_t &i_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_,const Bool_t &i_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_,const Bool_t &i_Photon90_R9Id90_HE10_IsoM_,const Bool_t &i_Photon120_R9Id90_HE10_IsoM_,const Bool_t &i_Photon165_R9Id90_HE10_IsoM_,const Bool_t &i_Photon90_CaloIdL_PFHT700_,const Bool_t &i_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_,const Bool_t &i_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_,const Bool_t &i_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_,const Bool_t &i_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_,const Bool_t &i_Photon35_TwoProngs35_,const Bool_t &i_IsoMu24_TwoProngs35_,const Bool_t &i_Dimuon0_Jpsi_L1_NoOS_,const Bool_t &i_Dimuon0_Jpsi_NoVertexing_NoOS_,const Bool_t &i_Dimuon0_Jpsi_,const Bool_t &i_Dimuon0_Jpsi_NoVertexing_,const Bool_t &i_Dimuon0_Jpsi_L1_4R_0er1p5R_,const Bool_t &i_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_,const Bool_t &i_Dimuon0_Jpsi3p5_Muon2_,const Bool_t &i_Dimuon0_Upsilon_L1_4p5_,const Bool_t &i_Dimuon0_Upsilon_L1_5_,const Bool_t &i_Dimuon0_Upsilon_L1_4p5NoOS_,const Bool_t &i_Dimuon0_Upsilon_L1_4p5er2p0_,const Bool_t &i_Dimuon0_Upsilon_L1_4p5er2p0M_,const Bool_t &i_Dimuon0_Upsilon_NoVertexing_,const Bool_t &i_Dimuon0_Upsilon_L1_5M_,const Bool_t &i_Dimuon0_LowMass_L1_0er1p5R_,const Bool_t &i_Dimuon0_LowMass_L1_0er1p5_,const Bool_t &i_Dimuon0_LowMass_,const Bool_t &i_Dimuon0_LowMass_L1_4_,const Bool_t &i_Dimuon0_LowMass_L1_4R_,const Bool_t &i_Dimuon0_LowMass_L1_TM530_,const Bool_t &i_Dimuon0_Upsilon_Muon_L1_TM0_,const Bool_t &i_Dimuon0_Upsilon_Muon_NoL1Mass_,const Bool_t &i_TripleMu_5_3_3_Mass3p8_DZ_,const Bool_t &i_TripleMu_10_5_5_DZ_,const Bool_t &i_TripleMu_12_10_5_,const Bool_t &i_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_,const Bool_t &i_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_,const Bool_t &i_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_,const Bool_t &i_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_,const Bool_t &i_DoubleMu3_DZ_PFMET50_PFMHT60_,const Bool_t &i_DoubleMu3_DZ_PFMET70_PFMHT70_,const Bool_t &i_DoubleMu3_DZ_PFMET90_PFMHT90_,const Bool_t &i_DoubleMu3_Trk_Tau3mu_NoL1Mass_,const Bool_t &i_DoubleMu4_Jpsi_Displaced_,const Bool_t &i_DoubleMu4_Jpsi_NoVertexing_,const Bool_t &i_DoubleMu4_JpsiTrkTrk_Displaced_,const Bool_t &i_DoubleMu43NoFiltersNoVtx_,const Bool_t &i_DoubleMu48NoFiltersNoVtx_,const Bool_t &i_Mu43NoFiltersNoVtx_Photon43_CaloIdL_,const Bool_t &i_Mu48NoFiltersNoVtx_Photon48_CaloIdL_,const Bool_t &i_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_,const Bool_t &i_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_,const Bool_t &i_DoubleMu33NoFiltersNoVtxDisplaced_,const Bool_t &i_DoubleMu40NoFiltersNoVtxDisplaced_,const Bool_t &i_DoubleMu20_7_Mass0to30_L1_DM4_,const Bool_t &i_DoubleMu20_7_Mass0to30_L1_DM4EG_,const Bool_t &i_HT425_,const Bool_t &i_HT430_DisplacedDijet40_DisplacedTrack_,const Bool_t &i_HT500_DisplacedDijet40_DisplacedTrack_,const Bool_t &i_HT430_DisplacedDijet60_DisplacedTrack_,const Bool_t &i_HT400_DisplacedDijet40_DisplacedTrack_,const Bool_t &i_HT650_DisplacedDijet60_Inclusive_,const Bool_t &i_HT550_DisplacedDijet60_Inclusive_,const Bool_t &i_DiJet110_35_Mjj650_PFMET110_,const Bool_t &i_DiJet110_35_Mjj650_PFMET120_,const Bool_t &i_DiJet110_35_Mjj650_PFMET130_,const Bool_t &i_TripleJet110_35_35_Mjj650_PFMET110_,const Bool_t &i_TripleJet110_35_35_Mjj650_PFMET120_,const Bool_t &i_TripleJet110_35_35_Mjj650_PFMET130_,const Bool_t &i_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_,const Bool_t &i_Ele28_eta2p1_WPTight_Gsf_HT150_,const Bool_t &i_Ele28_HighEta_SC20_Mass55_,const Bool_t &i_DoubleMu20_7_Mass0to30_Photon23_,const Bool_t &i_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_,const Bool_t &i_Ele15_IsoVVVL_PFHT450_PFMET50_,const Bool_t &i_Ele15_IsoVVVL_PFHT450_,const Bool_t &i_Ele50_IsoVVVL_PFHT450_,const Bool_t &i_Ele15_IsoVVVL_PFHT600_,const Bool_t &i_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_,const Bool_t &i_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_,const Bool_t &i_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_,const Bool_t &i_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_,const Bool_t &i_Mu15_IsoVVVL_PFHT450_PFMET50_,const Bool_t &i_Mu15_IsoVVVL_PFHT450_,const Bool_t &i_Mu50_IsoVVVL_PFHT450_,const Bool_t &i_Mu15_IsoVVVL_PFHT600_,const Bool_t &i_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_,const Bool_t &i_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_,const Bool_t &i_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_,const Bool_t &i_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_,const Bool_t &i_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_,const Bool_t &i_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_,const Bool_t &i_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_,const Bool_t &i_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_,const Bool_t &i_Dimuon10_PsiPrime_Barrel_Seagulls_,const Bool_t &i_Dimuon20_Jpsi_Barrel_Seagulls_,const Bool_t &i_Dimuon12_Upsilon_y1p4_,const Bool_t &i_Dimuon14_Phi_Barrel_Seagulls_,const Bool_t &i_Dimuon18_PsiPrime_,const Bool_t &i_Dimuon25_Jpsi_,const Bool_t &i_Dimuon18_PsiPrime_noCorrL1_,const Bool_t &i_Dimuon24_Upsilon_noCorrL1_,const Bool_t &i_Dimuon24_Phi_noCorrL1_,const Bool_t &i_Dimuon25_Jpsi_noCorrL1_,const Bool_t &i_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_,const Bool_t &i_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_,const Bool_t &i_DiMu9_Ele9_CaloIdL_TrackIdL_,const Bool_t &i_DoubleIsoMu20_eta2p1_,const Bool_t &i_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_,const Bool_t &i_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_,const Bool_t &i_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_,const Bool_t &i_Mu8_,const Bool_t &i_Mu17_,const Bool_t &i_Mu19_,const Bool_t &i_Mu17_Photon30_IsoCaloId_,const Bool_t &i_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_,const Bool_t &i_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_,const Bool_t &i_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_,const Bool_t &i_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_,const Bool_t &i_Ele8_CaloIdM_TrackIdM_PFJet30_,const Bool_t &i_Ele17_CaloIdM_TrackIdM_PFJet30_,const Bool_t &i_Ele23_CaloIdM_TrackIdM_PFJet30_,const Bool_t &i_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_,const Bool_t &i_Ele115_CaloIdVT_GsfTrkIdT_,const Bool_t &i_Ele135_CaloIdVT_GsfTrkIdT_,const Bool_t &i_Ele145_CaloIdVT_GsfTrkIdT_,const Bool_t &i_Ele200_CaloIdVT_GsfTrkIdT_,const Bool_t &i_Ele250_CaloIdVT_GsfTrkIdT_,const Bool_t &i_Ele300_CaloIdVT_GsfTrkIdT_,const Bool_t &i_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_,const Bool_t &i_PFHT330PT30_QuadPFJet_75_60_45_40_,const Bool_t &i_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_,const Bool_t &i_PFHT400_SixPFJet32_,const Bool_t &i_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_,const Bool_t &i_PFHT450_SixPFJet36_,const Bool_t &i_PFHT350_,const Bool_t &i_PFHT350MinPFJet15_,const Bool_t &i_Photon60_R9Id90_CaloIdL_IsoL_,const Bool_t &i_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_,const Bool_t &i_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_,const Bool_t &i_ECALHT800_,const Bool_t &i_DiSC30_18_EIso_AND_HE_Mass70_,const Bool_t &i_Physics_,const Bool_t &i_Physics_part0_,const Bool_t &i_Physics_part1_,const Bool_t &i_Physics_part2_,const Bool_t &i_Physics_part3_,const Bool_t &i_Physics_part4_,const Bool_t &i_Physics_part5_,const Bool_t &i_Physics_part6_,const Bool_t &i_Physics_part7_,const Bool_t &i_Random_,const Bool_t &i_ZeroBias_,const Bool_t &i_ZeroBias_Alignment_,const Bool_t &i_ZeroBias_part0_,const Bool_t &i_ZeroBias_part1_,const Bool_t &i_ZeroBias_part2_,const Bool_t &i_ZeroBias_part3_,const Bool_t &i_ZeroBias_part4_,const Bool_t &i_ZeroBias_part5_,const Bool_t &i_ZeroBias_part6_,const Bool_t &i_ZeroBias_part7_,const Bool_t &i_AK4CaloJet30_,const Bool_t &i_AK4CaloJet40_,const Bool_t &i_AK4CaloJet50_,const Bool_t &i_AK4CaloJet80_,const Bool_t &i_AK4CaloJet100_,const Bool_t &i_AK4CaloJet120_,const Bool_t &i_AK4PFJet30_,const Bool_t &i_AK4PFJet50_,const Bool_t &i_AK4PFJet80_,const Bool_t &i_AK4PFJet100_,const Bool_t &i_AK4PFJet120_,const Bool_t &i_SinglePhoton10_Eta3p1ForPPRef_,const Bool_t &i_SinglePhoton20_Eta3p1ForPPRef_,const Bool_t &i_SinglePhoton30_Eta3p1ForPPRef_,const Bool_t &i_Photon20_HoverELoose_,const Bool_t &i_Photon30_HoverELoose_,const Bool_t &i_EcalCalibration_,const Bool_t &i_HcalCalibration_,const Bool_t &i_L1UnpairedBunchBptxMinus_,const Bool_t &i_L1UnpairedBunchBptxPlus_,const Bool_t &i_L1NotBptxOR_,const Bool_t &i_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_,const Bool_t &i_CDC_L2cosmic_5_er1p0_,const Bool_t &i_CDC_L2cosmic_5p5_er1p0_,const Bool_t &i_HcalNZS_,const Bool_t &i_HcalPhiSym_,const Bool_t &i_HcalIsolatedbunch_,const Bool_t &i_IsoTrackHB_,const Bool_t &i_IsoTrackHE_,const Bool_t &i_ZeroBias_FirstCollisionAfterAbortGap_,const Bool_t &i_ZeroBias_IsolatedBunches_,const Bool_t &i_ZeroBias_FirstCollisionInTrain_,const Bool_t &i_ZeroBias_LastCollisionInTrain_,const Bool_t &i_ZeroBias_FirstBXAfterTrain_,const Bool_t &i_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_,const Bool_t &i_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_,const Bool_t &i_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_,const Bool_t &i_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_,const Bool_t &i_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_,const Bool_t &i_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_,const Bool_t &i_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_,const Bool_t &i_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_,const Bool_t &i_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_,const Bool_t &i_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_,const Bool_t &i_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_,const Bool_t &i_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_,const Bool_t &i_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_,const Bool_t &i_Rsq0p35_,const Bool_t &i_Rsq0p40_,const Bool_t &i_RsqMR300_Rsq0p09_MR200_,const Bool_t &i_RsqMR320_Rsq0p09_MR200_,const Bool_t &i_RsqMR300_Rsq0p09_MR200_4jet_,const Bool_t &i_RsqMR320_Rsq0p09_MR200_4jet_,const Bool_t &i_IsoMu27_MET90_,const Bool_t &i_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_,const Bool_t &i_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_,const Bool_t &i_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_,const Bool_t &i_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_,const Bool_t &i_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_,const Bool_t &i_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_,const Bool_t &i_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_,const Bool_t &i_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_,const Bool_t &i_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_,const Bool_t &i_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_,const Bool_t &i_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_,const Bool_t &i_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_,const Bool_t &i_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_,const Bool_t &i_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_,const Bool_t &i_PFMET100_PFMHT100_IDTight_PFHT60_,const Bool_t &i_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_,const Bool_t &i_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_,const Bool_t &i_Mu18_Mu9_SameSign_,const Bool_t &i_Mu18_Mu9_SameSign_DZ_,const Bool_t &i_Mu18_Mu9_,const Bool_t &i_Mu18_Mu9_DZ_,const Bool_t &i_Mu20_Mu10_SameSign_,const Bool_t &i_Mu20_Mu10_SameSign_DZ_,const Bool_t &i_Mu20_Mu10_,const Bool_t &i_Mu20_Mu10_DZ_,const Bool_t &i_Mu23_Mu12_SameSign_,const Bool_t &i_Mu23_Mu12_SameSign_DZ_,const Bool_t &i_Mu23_Mu12_,const Bool_t &i_Mu23_Mu12_DZ_,const Bool_t &i_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_,const Bool_t &i_DoubleMu2_Jpsi_DoubleTkMu0_Phi_,const Bool_t &i_DoubleMu3_DCA_PFMET50_PFMHT60_,const Bool_t &i_TripleMu_5_3_3_Mass3p8_DCA_,const Bool_t &i_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_,const Bool_t &i_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_,const Bool_t &i_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_,const Bool_t &i_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_,const Bool_t &i_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_,const Bool_t &i_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_,const Bool_t &i_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_,const Bool_t &i_QuadPFJet98_83_71_15_,const Bool_t &i_QuadPFJet103_88_75_15_,const Bool_t &i_QuadPFJet105_88_76_15_,const Bool_t &i_QuadPFJet111_90_80_15_,const Bool_t &i_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_,const Bool_t &i_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_,const Bool_t &i_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_,const Bool_t &i_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_,const Bool_t &i_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_,const Bool_t &i_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_,const Bool_t &i_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_,const Bool_t &i_Mu12_IP6_part0_,const Bool_t &i_Mu12_IP6_part1_,const Bool_t &i_Mu12_IP6_part2_,const Bool_t &i_Mu12_IP6_part3_,const Bool_t &i_Mu12_IP6_part4_,const Bool_t &i_Mu9_IP5_part0_,const Bool_t &i_Mu9_IP5_part1_,const Bool_t &i_Mu9_IP5_part2_,const Bool_t &i_Mu9_IP5_part3_,const Bool_t &i_Mu9_IP5_part4_,const Bool_t &i_Mu7_IP4_part0_,const Bool_t &i_Mu7_IP4_part1_,const Bool_t &i_Mu7_IP4_part2_,const Bool_t &i_Mu7_IP4_part3_,const Bool_t &i_Mu7_IP4_part4_,const Bool_t &i_Mu9_IP4_part0_,const Bool_t &i_Mu9_IP4_part1_,const Bool_t &i_Mu9_IP4_part2_,const Bool_t &i_Mu9_IP4_part3_,const Bool_t &i_Mu9_IP4_part4_,const Bool_t &i_Mu8_IP5_part0_,const Bool_t &i_Mu8_IP5_part1_,const Bool_t &i_Mu8_IP5_part2_,const Bool_t &i_Mu8_IP5_part3_,const Bool_t &i_Mu8_IP5_part4_,const Bool_t &i_Mu8_IP6_part0_,const Bool_t &i_Mu8_IP6_part1_,const Bool_t &i_Mu8_IP6_part2_,const Bool_t &i_Mu8_IP6_part3_,const Bool_t &i_Mu8_IP6_part4_,const Bool_t &i_Mu9_IP6_part0_,const Bool_t &i_Mu9_IP6_part1_,const Bool_t &i_Mu9_IP6_part2_,const Bool_t &i_Mu9_IP6_part3_,const Bool_t &i_Mu9_IP6_part4_,const Bool_t &i_Mu8_IP3_part0_,const Bool_t &i_Mu8_IP3_part1_,const Bool_t &i_Mu8_IP3_part2_,const Bool_t &i_Mu8_IP3_part3_,const Bool_t &i_Mu8_IP3_part4_,const Bool_t &i_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_,const Bool_t &i_TrkMu6NoFiltersNoVtx_,const Bool_t &i_TrkMu16NoFiltersNoVtx_,const Bool_t &i_DoubleTrkMu_16_6_NoFiltersNoVtx_,const Bool_t &i_HLTriggerFinalPath_):
//    
//  {}
  Hlt():
    HLTriggerFirstPath_(0),
    AK8PFJet360_TrimMass30_(0),
    AK8PFJet380_TrimMass30_(0),
    AK8PFJet400_TrimMass30_(0),
    AK8PFJet420_TrimMass30_(0),
    AK8PFHT750_TrimMass50_(0),
    AK8PFHT800_TrimMass50_(0),
    AK8PFHT850_TrimMass50_(0),
    AK8PFHT900_TrimMass50_(0),
    CaloJet500_NoJetID_(0),
    CaloJet550_NoJetID_(0),
    DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_(0),
    DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_(0),
    Trimuon5_3p5_2_Upsilon_Muon_(0),
    TrimuonOpen_5_3p5_2_Upsilon_Muon_(0),
    DoubleEle25_CaloIdL_MW_(0),
    DoubleEle27_CaloIdL_MW_(0),
    DoubleEle33_CaloIdL_MW_(0),
    DoubleEle24_eta2p1_WPTight_Gsf_(0),
    DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_(0),
    DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_(0),
    Ele27_Ele37_CaloIdL_MW_(0),
    Mu27_Ele37_CaloIdL_MW_(0),
    Mu37_Ele27_CaloIdL_MW_(0),
    Mu37_TkMu27_(0),
    DoubleMu4_3_Bs_(0),
    DoubleMu4_3_Jpsi_(0),
    DoubleMu4_JpsiTrk_Displaced_(0),
    DoubleMu4_LowMassNonResonantTrk_Displaced_(0),
    DoubleMu3_Trk_Tau3mu_(0),
    DoubleMu3_TkMu_DsTau3Mu_(0),
    DoubleMu4_PsiPrimeTrk_Displaced_(0),
    DoubleMu4_Mass3p8_DZ_PFHT350_(0),
    Mu3_PFJet40_(0),
    Mu7p5_L2Mu2_Jpsi_(0),
    Mu7p5_L2Mu2_Upsilon_(0),
    Mu7p5_Track2_Jpsi_(0),
    Mu7p5_Track3p5_Jpsi_(0),
    Mu7p5_Track7_Jpsi_(0),
    Mu7p5_Track2_Upsilon_(0),
    Mu7p5_Track3p5_Upsilon_(0),
    Mu7p5_Track7_Upsilon_(0),
    Mu3_L1SingleMu5orSingleMu7_(0),
    DoublePhoton33_CaloIdL_(0),
    DoublePhoton70_(0),
    DoublePhoton85_(0),
    Ele20_WPTight_Gsf_(0),
    Ele15_WPLoose_Gsf_(0),
    Ele17_WPLoose_Gsf_(0),
    Ele20_WPLoose_Gsf_(0),
    Ele20_eta2p1_WPLoose_Gsf_(0),
    DiEle27_WPTightCaloOnly_L1DoubleEG_(0),
    Ele27_WPTight_Gsf_(0),
    Ele28_WPTight_Gsf_(0),
    Ele30_WPTight_Gsf_(0),
    Ele32_WPTight_Gsf_(0),
    Ele35_WPTight_Gsf_(0),
    Ele35_WPTight_Gsf_L1EGMT_(0),
    Ele38_WPTight_Gsf_(0),
    Ele40_WPTight_Gsf_(0),
    Ele32_WPTight_Gsf_L1DoubleEG_(0),
    Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_(0),
    Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_(0),
    Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_(0),
    Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_(0),
    Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_(0),
    Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_(0),
    HT450_Beamspot_(0),
    HT300_Beamspot_(0),
    ZeroBias_Beamspot_(0),
    IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_(0),
    IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_(0),
    IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_(0),
    IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_(0),
    IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_(0),
    IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_(0),
    IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_(0),
    IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_(0),
    IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_(0),
    IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_(0),
    IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_(0),
    IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_(0),
    IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_(0),
    IsoMu20_(0),
    IsoMu24_(0),
    IsoMu24_eta2p1_(0),
    IsoMu27_(0),
    IsoMu30_(0),
    UncorrectedJetE30_NoBPTX_(0),
    UncorrectedJetE30_NoBPTX3BX_(0),
    UncorrectedJetE60_NoBPTX3BX_(0),
    UncorrectedJetE70_NoBPTX3BX_(0),
    L1SingleMu18_(0),
    L1SingleMu25_(0),
    L2Mu10_(0),
    L2Mu10_NoVertex_NoBPTX3BX_(0),
    L2Mu10_NoVertex_NoBPTX_(0),
    L2Mu45_NoVertex_3Sta_NoBPTX3BX_(0),
    L2Mu40_NoVertex_3Sta_NoBPTX3BX_(0),
    L2Mu50_(0),
    L2Mu23NoVtx_2Cha_(0),
    L2Mu23NoVtx_2Cha_CosmicSeed_(0),
    DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_(0),
    DoubleL2Mu30NoVtx_2Cha_Eta2p4_(0),
    DoubleL2Mu50_(0),
    DoubleL2Mu23NoVtx_2Cha_CosmicSeed_(0),
    DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_(0),
    DoubleL2Mu25NoVtx_2Cha_CosmicSeed_(0),
    DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_(0),
    DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_(0),
    DoubleL2Mu23NoVtx_2Cha_(0),
    DoubleL2Mu23NoVtx_2Cha_NoL2Matched_(0),
    DoubleL2Mu25NoVtx_2Cha_(0),
    DoubleL2Mu25NoVtx_2Cha_NoL2Matched_(0),
    DoubleL2Mu25NoVtx_2Cha_Eta2p4_(0),
    Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_(0),
    Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_(0),
    Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_(0),
    Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_(0),
    Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_(0),
    Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_(0),
    Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_(0),
    Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_(0),
    Mu25_TkMu0_Onia_(0),
    Mu30_TkMu0_Psi_(0),
    Mu30_TkMu0_Upsilon_(0),
    Mu20_TkMu0_Phi_(0),
    Mu25_TkMu0_Phi_(0),
    Mu12_(0),
    Mu15_(0),
    Mu20_(0),
    Mu27_(0),
    Mu50_(0),
    Mu55_(0),
    OldMu100_(0),
    TkMu100_(0),
    DiPFJetAve40_(0),
    DiPFJetAve60_(0),
    DiPFJetAve80_(0),
    DiPFJetAve140_(0),
    DiPFJetAve200_(0),
    DiPFJetAve260_(0),
    DiPFJetAve320_(0),
    DiPFJetAve400_(0),
    DiPFJetAve500_(0),
    DiPFJetAve60_HFJEC_(0),
    DiPFJetAve80_HFJEC_(0),
    DiPFJetAve100_HFJEC_(0),
    DiPFJetAve160_HFJEC_(0),
    DiPFJetAve220_HFJEC_(0),
    DiPFJetAve300_HFJEC_(0),
    AK8PFJet15_(0),
    AK8PFJet25_(0),
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
    AK8PFJet550_(0),
    PFJet15_(0),
    PFJet25_(0),
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
    PFJet550_(0),
    PFJetFwd15_(0),
    PFJetFwd25_(0),
    PFJetFwd40_(0),
    PFJetFwd60_(0),
    PFJetFwd80_(0),
    PFJetFwd140_(0),
    PFJetFwd200_(0),
    PFJetFwd260_(0),
    PFJetFwd320_(0),
    PFJetFwd400_(0),
    PFJetFwd450_(0),
    PFJetFwd500_(0),
    AK8PFJetFwd15_(0),
    AK8PFJetFwd25_(0),
    AK8PFJetFwd40_(0),
    AK8PFJetFwd60_(0),
    AK8PFJetFwd80_(0),
    AK8PFJetFwd140_(0),
    AK8PFJetFwd200_(0),
    AK8PFJetFwd260_(0),
    AK8PFJetFwd320_(0),
    AK8PFJetFwd400_(0),
    AK8PFJetFwd450_(0),
    AK8PFJetFwd500_(0),
    PFHT180_(0),
    PFHT250_(0),
    PFHT370_(0),
    PFHT430_(0),
    PFHT510_(0),
    PFHT590_(0),
    PFHT680_(0),
    PFHT780_(0),
    PFHT890_(0),
    PFHT1050_(0),
    PFHT500_PFMET100_PFMHT100_IDTight_(0),
    PFHT500_PFMET110_PFMHT110_IDTight_(0),
    PFHT700_PFMET85_PFMHT85_IDTight_(0),
    PFHT700_PFMET95_PFMHT95_IDTight_(0),
    PFHT800_PFMET75_PFMHT75_IDTight_(0),
    PFHT800_PFMET85_PFMHT85_IDTight_(0),
    PFMET110_PFMHT110_IDTight_(0),
    PFMET120_PFMHT120_IDTight_(0),
    PFMET130_PFMHT130_IDTight_(0),
    PFMET140_PFMHT140_IDTight_(0),
    PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_(0),
    PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_(0),
    PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_(0),
    PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_(0),
    PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_(0),
    PFMET120_PFMHT120_IDTight_PFHT60_(0),
    PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_(0),
    PFMETTypeOne120_PFMHT120_IDTight_PFHT60_(0),
    PFMETTypeOne110_PFMHT110_IDTight_(0),
    PFMETTypeOne120_PFMHT120_IDTight_(0),
    PFMETTypeOne130_PFMHT130_IDTight_(0),
    PFMETTypeOne140_PFMHT140_IDTight_(0),
    PFMETNoMu110_PFMHTNoMu110_IDTight_(0),
    PFMETNoMu120_PFMHTNoMu120_IDTight_(0),
    PFMETNoMu130_PFMHTNoMu130_IDTight_(0),
    PFMETNoMu140_PFMHTNoMu140_IDTight_(0),
    MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_(0),
    MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_(0),
    MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_(0),
    MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_(0),
    L1ETMHadSeeds_(0),
    CaloMHT90_(0),
    CaloMET80_NotCleaned_(0),
    CaloMET90_NotCleaned_(0),
    CaloMET100_NotCleaned_(0),
    CaloMET110_NotCleaned_(0),
    CaloMET250_NotCleaned_(0),
    CaloMET70_HBHECleaned_(0),
    CaloMET80_HBHECleaned_(0),
    CaloMET90_HBHECleaned_(0),
    CaloMET100_HBHECleaned_(0),
    CaloMET250_HBHECleaned_(0),
    CaloMET300_HBHECleaned_(0),
    CaloMET350_HBHECleaned_(0),
    PFMET200_NotCleaned_(0),
    PFMET200_HBHECleaned_(0),
    PFMET250_HBHECleaned_(0),
    PFMET300_HBHECleaned_(0),
    PFMET200_HBHE_BeamHaloCleaned_(0),
    PFMETTypeOne200_HBHE_BeamHaloCleaned_(0),
    MET105_IsoTrk50_(0),
    MET120_IsoTrk50_(0),
    SingleJet30_Mu12_SinglePFJet40_(0),
    Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_(0),
    Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_(0),
    Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_(0),
    Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_(0),
    Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    DoublePFJets40_CaloBTagDeepCSV_p71_(0),
    DoublePFJets100_CaloBTagDeepCSV_p71_(0),
    DoublePFJets200_CaloBTagDeepCSV_p71_(0),
    DoublePFJets350_CaloBTagDeepCSV_p71_(0),
    DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    Photon300_NoHE_(0),
    Mu8_TrkIsoVVL_(0),
    Mu8_DiEle12_CaloIdL_TrackIdL_DZ_(0),
    Mu8_DiEle12_CaloIdL_TrackIdL_(0),
    Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_(0),
    Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_(0),
    Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_(0),
    Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_(0),
    Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_(0),
    Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_(0),
    Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(0),
    Mu17_TrkIsoVVL_(0),
    Mu19_TrkIsoVVL_(0),
    BTagMu_AK4DiJet20_Mu5_(0),
    BTagMu_AK4DiJet40_Mu5_(0),
    BTagMu_AK4DiJet70_Mu5_(0),
    BTagMu_AK4DiJet110_Mu5_(0),
    BTagMu_AK4DiJet170_Mu5_(0),
    BTagMu_AK4Jet300_Mu5_(0),
    BTagMu_AK8DiJet170_Mu5_(0),
    BTagMu_AK8Jet170_DoubleMu5_(0),
    BTagMu_AK8Jet300_Mu5_(0),
    BTagMu_AK4DiJet20_Mu5_noalgo_(0),
    BTagMu_AK4DiJet40_Mu5_noalgo_(0),
    BTagMu_AK4DiJet70_Mu5_noalgo_(0),
    BTagMu_AK4DiJet110_Mu5_noalgo_(0),
    BTagMu_AK4DiJet170_Mu5_noalgo_(0),
    BTagMu_AK4Jet300_Mu5_noalgo_(0),
    BTagMu_AK8DiJet170_Mu5_noalgo_(0),
    BTagMu_AK8Jet170_DoubleMu5_noalgo_(0),
    BTagMu_AK8Jet300_Mu5_noalgo_(0),
    Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_(0),
    Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(0),
    Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    Mu12_DoublePhoton20_(0),
    TriplePhoton_20_20_20_CaloIdLV2_(0),
    TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_(0),
    TriplePhoton_30_30_10_CaloIdLV2_(0),
    TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_(0),
    TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_(0),
    Photon20_(0),
    Photon33_(0),
    Photon50_(0),
    Photon75_(0),
    Photon90_(0),
    Photon120_(0),
    Photon150_(0),
    Photon175_(0),
    Photon200_(0),
    Photon100EB_TightID_TightIso_(0),
    Photon110EB_TightID_TightIso_(0),
    Photon120EB_TightID_TightIso_(0),
    Photon100EBHE10_(0),
    Photon100EEHE10_(0),
    Photon100EE_TightID_TightIso_(0),
    Photon50_R9Id90_HE10_IsoM_(0),
    Photon75_R9Id90_HE10_IsoM_(0),
    Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_(0),
    Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_(0),
    Photon90_R9Id90_HE10_IsoM_(0),
    Photon120_R9Id90_HE10_IsoM_(0),
    Photon165_R9Id90_HE10_IsoM_(0),
    Photon90_CaloIdL_PFHT700_(0),
    Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_(0),
    Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_(0),
    Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_(0),
    Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_(0),
    Photon35_TwoProngs35_(0),
    IsoMu24_TwoProngs35_(0),
    Dimuon0_Jpsi_L1_NoOS_(0),
    Dimuon0_Jpsi_NoVertexing_NoOS_(0),
    Dimuon0_Jpsi_(0),
    Dimuon0_Jpsi_NoVertexing_(0),
    Dimuon0_Jpsi_L1_4R_0er1p5R_(0),
    Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_(0),
    Dimuon0_Jpsi3p5_Muon2_(0),
    Dimuon0_Upsilon_L1_4p5_(0),
    Dimuon0_Upsilon_L1_5_(0),
    Dimuon0_Upsilon_L1_4p5NoOS_(0),
    Dimuon0_Upsilon_L1_4p5er2p0_(0),
    Dimuon0_Upsilon_L1_4p5er2p0M_(0),
    Dimuon0_Upsilon_NoVertexing_(0),
    Dimuon0_Upsilon_L1_5M_(0),
    Dimuon0_LowMass_L1_0er1p5R_(0),
    Dimuon0_LowMass_L1_0er1p5_(0),
    Dimuon0_LowMass_(0),
    Dimuon0_LowMass_L1_4_(0),
    Dimuon0_LowMass_L1_4R_(0),
    Dimuon0_LowMass_L1_TM530_(0),
    Dimuon0_Upsilon_Muon_L1_TM0_(0),
    Dimuon0_Upsilon_Muon_NoL1Mass_(0),
    TripleMu_5_3_3_Mass3p8_DZ_(0),
    TripleMu_10_5_5_DZ_(0),
    TripleMu_12_10_5_(0),
    Tau3Mu_Mu7_Mu1_TkMu1_Tau15_(0),
    Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_(0),
    Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_(0),
    Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_(0),
    DoubleMu3_DZ_PFMET50_PFMHT60_(0),
    DoubleMu3_DZ_PFMET70_PFMHT70_(0),
    DoubleMu3_DZ_PFMET90_PFMHT90_(0),
    DoubleMu3_Trk_Tau3mu_NoL1Mass_(0),
    DoubleMu4_Jpsi_Displaced_(0),
    DoubleMu4_Jpsi_NoVertexing_(0),
    DoubleMu4_JpsiTrkTrk_Displaced_(0),
    DoubleMu43NoFiltersNoVtx_(0),
    DoubleMu48NoFiltersNoVtx_(0),
    Mu43NoFiltersNoVtx_Photon43_CaloIdL_(0),
    Mu48NoFiltersNoVtx_Photon48_CaloIdL_(0),
    Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_(0),
    Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_(0),
    DoubleMu33NoFiltersNoVtxDisplaced_(0),
    DoubleMu40NoFiltersNoVtxDisplaced_(0),
    DoubleMu20_7_Mass0to30_L1_DM4_(0),
    DoubleMu20_7_Mass0to30_L1_DM4EG_(0),
    HT425_(0),
    HT430_DisplacedDijet40_DisplacedTrack_(0),
    HT500_DisplacedDijet40_DisplacedTrack_(0),
    HT430_DisplacedDijet60_DisplacedTrack_(0),
    HT400_DisplacedDijet40_DisplacedTrack_(0),
    HT650_DisplacedDijet60_Inclusive_(0),
    HT550_DisplacedDijet60_Inclusive_(0),
    DiJet110_35_Mjj650_PFMET110_(0),
    DiJet110_35_Mjj650_PFMET120_(0),
    DiJet110_35_Mjj650_PFMET130_(0),
    TripleJet110_35_35_Mjj650_PFMET110_(0),
    TripleJet110_35_35_Mjj650_PFMET120_(0),
    TripleJet110_35_35_Mjj650_PFMET130_(0),
    Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_(0),
    Ele28_eta2p1_WPTight_Gsf_HT150_(0),
    Ele28_HighEta_SC20_Mass55_(0),
    DoubleMu20_7_Mass0to30_Photon23_(0),
    Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_(0),
    Ele15_IsoVVVL_PFHT450_PFMET50_(0),
    Ele15_IsoVVVL_PFHT450_(0),
    Ele50_IsoVVVL_PFHT450_(0),
    Ele15_IsoVVVL_PFHT600_(0),
    Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_(0),
    Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_(0),
    Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_(0),
    Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_(0),
    Mu15_IsoVVVL_PFHT450_PFMET50_(0),
    Mu15_IsoVVVL_PFHT450_(0),
    Mu50_IsoVVVL_PFHT450_(0),
    Mu15_IsoVVVL_PFHT600_(0),
    Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_(0),
    Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_(0),
    Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_(0),
    Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_(0),
    Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_(0),
    Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_(0),
    Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_(0),
    Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_(0),
    Dimuon10_PsiPrime_Barrel_Seagulls_(0),
    Dimuon20_Jpsi_Barrel_Seagulls_(0),
    Dimuon12_Upsilon_y1p4_(0),
    Dimuon14_Phi_Barrel_Seagulls_(0),
    Dimuon18_PsiPrime_(0),
    Dimuon25_Jpsi_(0),
    Dimuon18_PsiPrime_noCorrL1_(0),
    Dimuon24_Upsilon_noCorrL1_(0),
    Dimuon24_Phi_noCorrL1_(0),
    Dimuon25_Jpsi_noCorrL1_(0),
    DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_(0),
    DiMu9_Ele9_CaloIdL_TrackIdL_DZ_(0),
    DiMu9_Ele9_CaloIdL_TrackIdL_(0),
    DoubleIsoMu20_eta2p1_(0),
    TrkMu12_DoubleTrkMu5NoFiltersNoVtx_(0),
    TrkMu16_DoubleTrkMu6NoFiltersNoVtx_(0),
    TrkMu17_DoubleTrkMu8NoFiltersNoVtx_(0),
    Mu8_(0),
    Mu17_(0),
    Mu19_(0),
    Mu17_Photon30_IsoCaloId_(0),
    Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    Ele8_CaloIdM_TrackIdM_PFJet30_(0),
    Ele17_CaloIdM_TrackIdM_PFJet30_(0),
    Ele23_CaloIdM_TrackIdM_PFJet30_(0),
    Ele50_CaloIdVT_GsfTrkIdT_PFJet165_(0),
    Ele115_CaloIdVT_GsfTrkIdT_(0),
    Ele135_CaloIdVT_GsfTrkIdT_(0),
    Ele145_CaloIdVT_GsfTrkIdT_(0),
    Ele200_CaloIdVT_GsfTrkIdT_(0),
    Ele250_CaloIdVT_GsfTrkIdT_(0),
    Ele300_CaloIdVT_GsfTrkIdT_(0),
    PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_(0),
    PFHT330PT30_QuadPFJet_75_60_45_40_(0),
    PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_(0),
    PFHT400_SixPFJet32_(0),
    PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_(0),
    PFHT450_SixPFJet36_(0),
    PFHT350_(0),
    PFHT350MinPFJet15_(0),
    Photon60_R9Id90_CaloIdL_IsoL_(0),
    Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_(0),
    Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_(0),
    ECALHT800_(0),
    DiSC30_18_EIso_AND_HE_Mass70_(0),
    Physics_(0),
    Physics_part0_(0),
    Physics_part1_(0),
    Physics_part2_(0),
    Physics_part3_(0),
    Physics_part4_(0),
    Physics_part5_(0),
    Physics_part6_(0),
    Physics_part7_(0),
    Random_(0),
    ZeroBias_(0),
    ZeroBias_Alignment_(0),
    ZeroBias_part0_(0),
    ZeroBias_part1_(0),
    ZeroBias_part2_(0),
    ZeroBias_part3_(0),
    ZeroBias_part4_(0),
    ZeroBias_part5_(0),
    ZeroBias_part6_(0),
    ZeroBias_part7_(0),
    AK4CaloJet30_(0),
    AK4CaloJet40_(0),
    AK4CaloJet50_(0),
    AK4CaloJet80_(0),
    AK4CaloJet100_(0),
    AK4CaloJet120_(0),
    AK4PFJet30_(0),
    AK4PFJet50_(0),
    AK4PFJet80_(0),
    AK4PFJet100_(0),
    AK4PFJet120_(0),
    SinglePhoton10_Eta3p1ForPPRef_(0),
    SinglePhoton20_Eta3p1ForPPRef_(0),
    SinglePhoton30_Eta3p1ForPPRef_(0),
    Photon20_HoverELoose_(0),
    Photon30_HoverELoose_(0),
    EcalCalibration_(0),
    HcalCalibration_(0),
    L1UnpairedBunchBptxMinus_(0),
    L1UnpairedBunchBptxPlus_(0),
    L1NotBptxOR_(0),
    L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_(0),
    CDC_L2cosmic_5_er1p0_(0),
    CDC_L2cosmic_5p5_er1p0_(0),
    HcalNZS_(0),
    HcalPhiSym_(0),
    HcalIsolatedbunch_(0),
    IsoTrackHB_(0),
    IsoTrackHE_(0),
    ZeroBias_FirstCollisionAfterAbortGap_(0),
    ZeroBias_IsolatedBunches_(0),
    ZeroBias_FirstCollisionInTrain_(0),
    ZeroBias_LastCollisionInTrain_(0),
    ZeroBias_FirstBXAfterTrain_(0),
    IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_(0),
    MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_(0),
    MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_(0),
    MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_(0),
    MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_(0),
    MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_(0),
    MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_(0),
    MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_(0),
    MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_(0),
    MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_(0),
    MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_(0),
    MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_(0),
    Ele16_Ele12_Ele8_CaloIdL_TrackIdL_(0),
    Rsq0p35_(0),
    Rsq0p40_(0),
    RsqMR300_Rsq0p09_MR200_(0),
    RsqMR320_Rsq0p09_MR200_(0),
    RsqMR300_Rsq0p09_MR200_4jet_(0),
    RsqMR320_Rsq0p09_MR200_4jet_(0),
    IsoMu27_MET90_(0),
    DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_(0),
    DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_(0),
    DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_(0),
    DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_(0),
    DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_(0),
    DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_(0),
    DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_(0),
    DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_(0),
    VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_(0),
    VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_(0),
    VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_(0),
    Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_(0),
    Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_(0),
    Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_(0),
    PFMET100_PFMHT100_IDTight_PFHT60_(0),
    PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_(0),
    PFMETTypeOne100_PFMHT100_IDTight_PFHT60_(0),
    Mu18_Mu9_SameSign_(0),
    Mu18_Mu9_SameSign_DZ_(0),
    Mu18_Mu9_(0),
    Mu18_Mu9_DZ_(0),
    Mu20_Mu10_SameSign_(0),
    Mu20_Mu10_SameSign_DZ_(0),
    Mu20_Mu10_(0),
    Mu20_Mu10_DZ_(0),
    Mu23_Mu12_SameSign_(0),
    Mu23_Mu12_SameSign_DZ_(0),
    Mu23_Mu12_(0),
    Mu23_Mu12_DZ_(0),
    DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_(0),
    DoubleMu2_Jpsi_DoubleTkMu0_Phi_(0),
    DoubleMu3_DCA_PFMET50_PFMHT60_(0),
    TripleMu_5_3_3_Mass3p8_DCA_(0),
    QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_(0),
    QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_(0),
    QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_(0),
    QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_(0),
    QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_(0),
    QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_(0),
    QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_(0),
    QuadPFJet98_83_71_15_(0),
    QuadPFJet103_88_75_15_(0),
    QuadPFJet105_88_76_15_(0),
    QuadPFJet111_90_80_15_(0),
    AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_(0),
    AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_(0),
    AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_(0),
    AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_(0),
    AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_(0),
    Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_(0),
    Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_(0),
    Mu12_IP6_part0_(0),
    Mu12_IP6_part1_(0),
    Mu12_IP6_part2_(0),
    Mu12_IP6_part3_(0),
    Mu12_IP6_part4_(0),
    Mu9_IP5_part0_(0),
    Mu9_IP5_part1_(0),
    Mu9_IP5_part2_(0),
    Mu9_IP5_part3_(0),
    Mu9_IP5_part4_(0),
    Mu7_IP4_part0_(0),
    Mu7_IP4_part1_(0),
    Mu7_IP4_part2_(0),
    Mu7_IP4_part3_(0),
    Mu7_IP4_part4_(0),
    Mu9_IP4_part0_(0),
    Mu9_IP4_part1_(0),
    Mu9_IP4_part2_(0),
    Mu9_IP4_part3_(0),
    Mu9_IP4_part4_(0),
    Mu8_IP5_part0_(0),
    Mu8_IP5_part1_(0),
    Mu8_IP5_part2_(0),
    Mu8_IP5_part3_(0),
    Mu8_IP5_part4_(0),
    Mu8_IP6_part0_(0),
    Mu8_IP6_part1_(0),
    Mu8_IP6_part2_(0),
    Mu8_IP6_part3_(0),
    Mu8_IP6_part4_(0),
    Mu9_IP6_part0_(0),
    Mu9_IP6_part1_(0),
    Mu9_IP6_part2_(0),
    Mu9_IP6_part3_(0),
    Mu9_IP6_part4_(0),
    Mu8_IP3_part0_(0),
    Mu8_IP3_part1_(0),
    Mu8_IP3_part2_(0),
    Mu8_IP3_part3_(0),
    Mu8_IP3_part4_(0),
    QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_(0),
    TrkMu6NoFiltersNoVtx_(0),
    TrkMu16NoFiltersNoVtx_(0),
    DoubleTrkMu_16_6_NoFiltersNoVtx_(0),
    HLTriggerFinalPath_(0)
  {}
  Bool_t HLTriggerFirstPath() const {return HLTriggerFirstPath_;}
  Bool_t AK8PFJet360_TrimMass30() const {return AK8PFJet360_TrimMass30_;}
  Bool_t AK8PFJet380_TrimMass30() const {return AK8PFJet380_TrimMass30_;}
  Bool_t AK8PFJet400_TrimMass30() const {return AK8PFJet400_TrimMass30_;}
  Bool_t AK8PFJet420_TrimMass30() const {return AK8PFJet420_TrimMass30_;}
  Bool_t AK8PFHT750_TrimMass50() const {return AK8PFHT750_TrimMass50_;}
  Bool_t AK8PFHT800_TrimMass50() const {return AK8PFHT800_TrimMass50_;}
  Bool_t AK8PFHT850_TrimMass50() const {return AK8PFHT850_TrimMass50_;}
  Bool_t AK8PFHT900_TrimMass50() const {return AK8PFHT900_TrimMass50_;}
  Bool_t CaloJet500_NoJetID() const {return CaloJet500_NoJetID_;}
  Bool_t CaloJet550_NoJetID() const {return CaloJet550_NoJetID_;}
  Bool_t DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL() const {return DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_;}
  Bool_t DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon() const {return DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_;}
  Bool_t Trimuon5_3p5_2_Upsilon_Muon() const {return Trimuon5_3p5_2_Upsilon_Muon_;}
  Bool_t TrimuonOpen_5_3p5_2_Upsilon_Muon() const {return TrimuonOpen_5_3p5_2_Upsilon_Muon_;}
  Bool_t DoubleEle25_CaloIdL_MW() const {return DoubleEle25_CaloIdL_MW_;}
  Bool_t DoubleEle27_CaloIdL_MW() const {return DoubleEle27_CaloIdL_MW_;}
  Bool_t DoubleEle33_CaloIdL_MW() const {return DoubleEle33_CaloIdL_MW_;}
  Bool_t DoubleEle24_eta2p1_WPTight_Gsf() const {return DoubleEle24_eta2p1_WPTight_Gsf_;}
  Bool_t DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350() const {return DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_;}
  Bool_t DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350() const {return DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_;}
  Bool_t Ele27_Ele37_CaloIdL_MW() const {return Ele27_Ele37_CaloIdL_MW_;}
  Bool_t Mu27_Ele37_CaloIdL_MW() const {return Mu27_Ele37_CaloIdL_MW_;}
  Bool_t Mu37_Ele27_CaloIdL_MW() const {return Mu37_Ele27_CaloIdL_MW_;}
  Bool_t Mu37_TkMu27() const {return Mu37_TkMu27_;}
  Bool_t DoubleMu4_3_Bs() const {return DoubleMu4_3_Bs_;}
  Bool_t DoubleMu4_3_Jpsi() const {return DoubleMu4_3_Jpsi_;}
  Bool_t DoubleMu4_JpsiTrk_Displaced() const {return DoubleMu4_JpsiTrk_Displaced_;}
  Bool_t DoubleMu4_LowMassNonResonantTrk_Displaced() const {return DoubleMu4_LowMassNonResonantTrk_Displaced_;}
  Bool_t DoubleMu3_Trk_Tau3mu() const {return DoubleMu3_Trk_Tau3mu_;}
  Bool_t DoubleMu3_TkMu_DsTau3Mu() const {return DoubleMu3_TkMu_DsTau3Mu_;}
  Bool_t DoubleMu4_PsiPrimeTrk_Displaced() const {return DoubleMu4_PsiPrimeTrk_Displaced_;}
  Bool_t DoubleMu4_Mass3p8_DZ_PFHT350() const {return DoubleMu4_Mass3p8_DZ_PFHT350_;}
  Bool_t Mu3_PFJet40() const {return Mu3_PFJet40_;}
  Bool_t Mu7p5_L2Mu2_Jpsi() const {return Mu7p5_L2Mu2_Jpsi_;}
  Bool_t Mu7p5_L2Mu2_Upsilon() const {return Mu7p5_L2Mu2_Upsilon_;}
  Bool_t Mu7p5_Track2_Jpsi() const {return Mu7p5_Track2_Jpsi_;}
  Bool_t Mu7p5_Track3p5_Jpsi() const {return Mu7p5_Track3p5_Jpsi_;}
  Bool_t Mu7p5_Track7_Jpsi() const {return Mu7p5_Track7_Jpsi_;}
  Bool_t Mu7p5_Track2_Upsilon() const {return Mu7p5_Track2_Upsilon_;}
  Bool_t Mu7p5_Track3p5_Upsilon() const {return Mu7p5_Track3p5_Upsilon_;}
  Bool_t Mu7p5_Track7_Upsilon() const {return Mu7p5_Track7_Upsilon_;}
  Bool_t Mu3_L1SingleMu5orSingleMu7() const {return Mu3_L1SingleMu5orSingleMu7_;}
  Bool_t DoublePhoton33_CaloIdL() const {return DoublePhoton33_CaloIdL_;}
  Bool_t DoublePhoton70() const {return DoublePhoton70_;}
  Bool_t DoublePhoton85() const {return DoublePhoton85_;}
  Bool_t Ele20_WPTight_Gsf() const {return Ele20_WPTight_Gsf_;}
  Bool_t Ele15_WPLoose_Gsf() const {return Ele15_WPLoose_Gsf_;}
  Bool_t Ele17_WPLoose_Gsf() const {return Ele17_WPLoose_Gsf_;}
  Bool_t Ele20_WPLoose_Gsf() const {return Ele20_WPLoose_Gsf_;}
  Bool_t Ele20_eta2p1_WPLoose_Gsf() const {return Ele20_eta2p1_WPLoose_Gsf_;}
  Bool_t DiEle27_WPTightCaloOnly_L1DoubleEG() const {return DiEle27_WPTightCaloOnly_L1DoubleEG_;}
  Bool_t Ele27_WPTight_Gsf() const {return Ele27_WPTight_Gsf_;}
  Bool_t Ele28_WPTight_Gsf() const {return Ele28_WPTight_Gsf_;}
  Bool_t Ele30_WPTight_Gsf() const {return Ele30_WPTight_Gsf_;}
  Bool_t Ele32_WPTight_Gsf() const {return Ele32_WPTight_Gsf_;}
  Bool_t Ele35_WPTight_Gsf() const {return Ele35_WPTight_Gsf_;}
  Bool_t Ele35_WPTight_Gsf_L1EGMT() const {return Ele35_WPTight_Gsf_L1EGMT_;}
  Bool_t Ele38_WPTight_Gsf() const {return Ele38_WPTight_Gsf_;}
  Bool_t Ele40_WPTight_Gsf() const {return Ele40_WPTight_Gsf_;}
  Bool_t Ele32_WPTight_Gsf_L1DoubleEG() const {return Ele32_WPTight_Gsf_L1DoubleEG_;}
  Bool_t Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1() const {return Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_;}
  Bool_t Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1() const {return Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_;}
  Bool_t Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1() const {return Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_;}
  Bool_t Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1() const {return Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;}
  Bool_t Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1() const {return Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;}
  Bool_t Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1() const {return Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;}
  Bool_t HT450_Beamspot() const {return HT450_Beamspot_;}
  Bool_t HT300_Beamspot() const {return HT300_Beamspot_;}
  Bool_t ZeroBias_Beamspot() const {return ZeroBias_Beamspot_;}
  Bool_t IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1() const {return IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_;}
  Bool_t IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1() const {return IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_;}
  Bool_t IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1() const {return IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_;}
  Bool_t IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1() const {return IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;}
  Bool_t IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1() const {return IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;}
  Bool_t IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1() const {return IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;}
  Bool_t IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1() const {return IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_;}
  Bool_t IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1() const {return IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_;}
  Bool_t IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1() const {return IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_;}
  Bool_t IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1() const {return IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_;}
  Bool_t IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1() const {return IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;}
  Bool_t IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1() const {return IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;}
  Bool_t IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1() const {return IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;}
  Bool_t IsoMu20() const {return IsoMu20_;}
  Bool_t IsoMu24() const {return IsoMu24_;}
  Bool_t IsoMu24_eta2p1() const {return IsoMu24_eta2p1_;}
  Bool_t IsoMu27() const {return IsoMu27_;}
  Bool_t IsoMu30() const {return IsoMu30_;}
  Bool_t UncorrectedJetE30_NoBPTX() const {return UncorrectedJetE30_NoBPTX_;}
  Bool_t UncorrectedJetE30_NoBPTX3BX() const {return UncorrectedJetE30_NoBPTX3BX_;}
  Bool_t UncorrectedJetE60_NoBPTX3BX() const {return UncorrectedJetE60_NoBPTX3BX_;}
  Bool_t UncorrectedJetE70_NoBPTX3BX() const {return UncorrectedJetE70_NoBPTX3BX_;}
  Bool_t L1SingleMu18() const {return L1SingleMu18_;}
  Bool_t L1SingleMu25() const {return L1SingleMu25_;}
  Bool_t L2Mu10() const {return L2Mu10_;}
  Bool_t L2Mu10_NoVertex_NoBPTX3BX() const {return L2Mu10_NoVertex_NoBPTX3BX_;}
  Bool_t L2Mu10_NoVertex_NoBPTX() const {return L2Mu10_NoVertex_NoBPTX_;}
  Bool_t L2Mu45_NoVertex_3Sta_NoBPTX3BX() const {return L2Mu45_NoVertex_3Sta_NoBPTX3BX_;}
  Bool_t L2Mu40_NoVertex_3Sta_NoBPTX3BX() const {return L2Mu40_NoVertex_3Sta_NoBPTX3BX_;}
  Bool_t L2Mu50() const {return L2Mu50_;}
  Bool_t L2Mu23NoVtx_2Cha() const {return L2Mu23NoVtx_2Cha_;}
  Bool_t L2Mu23NoVtx_2Cha_CosmicSeed() const {return L2Mu23NoVtx_2Cha_CosmicSeed_;}
  Bool_t DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4() const {return DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_;}
  Bool_t DoubleL2Mu30NoVtx_2Cha_Eta2p4() const {return DoubleL2Mu30NoVtx_2Cha_Eta2p4_;}
  Bool_t DoubleL2Mu50() const {return DoubleL2Mu50_;}
  Bool_t DoubleL2Mu23NoVtx_2Cha_CosmicSeed() const {return DoubleL2Mu23NoVtx_2Cha_CosmicSeed_;}
  Bool_t DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched() const {return DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_;}
  Bool_t DoubleL2Mu25NoVtx_2Cha_CosmicSeed() const {return DoubleL2Mu25NoVtx_2Cha_CosmicSeed_;}
  Bool_t DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched() const {return DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_;}
  Bool_t DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4() const {return DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_;}
  Bool_t DoubleL2Mu23NoVtx_2Cha() const {return DoubleL2Mu23NoVtx_2Cha_;}
  Bool_t DoubleL2Mu23NoVtx_2Cha_NoL2Matched() const {return DoubleL2Mu23NoVtx_2Cha_NoL2Matched_;}
  Bool_t DoubleL2Mu25NoVtx_2Cha() const {return DoubleL2Mu25NoVtx_2Cha_;}
  Bool_t DoubleL2Mu25NoVtx_2Cha_NoL2Matched() const {return DoubleL2Mu25NoVtx_2Cha_NoL2Matched_;}
  Bool_t DoubleL2Mu25NoVtx_2Cha_Eta2p4() const {return DoubleL2Mu25NoVtx_2Cha_Eta2p4_;}
  Bool_t Mu17_TrkIsoVVL_Mu8_TrkIsoVVL() const {return Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_;}
  Bool_t Mu19_TrkIsoVVL_Mu9_TrkIsoVVL() const {return Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_;}
  Bool_t Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() const {return Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_;}
  Bool_t Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ() const {return Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_;}
  Bool_t Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8() const {return Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_;}
  Bool_t Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8() const {return Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_;}
  Bool_t Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() const {return Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_;}
  Bool_t Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8() const {return Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_;}
  Bool_t Mu25_TkMu0_Onia() const {return Mu25_TkMu0_Onia_;}
  Bool_t Mu30_TkMu0_Psi() const {return Mu30_TkMu0_Psi_;}
  Bool_t Mu30_TkMu0_Upsilon() const {return Mu30_TkMu0_Upsilon_;}
  Bool_t Mu20_TkMu0_Phi() const {return Mu20_TkMu0_Phi_;}
  Bool_t Mu25_TkMu0_Phi() const {return Mu25_TkMu0_Phi_;}
  Bool_t Mu12() const {return Mu12_;}
  Bool_t Mu15() const {return Mu15_;}
  Bool_t Mu20() const {return Mu20_;}
  Bool_t Mu27() const {return Mu27_;}
  Bool_t Mu50() const {return Mu50_;}
  Bool_t Mu55() const {return Mu55_;}
  Bool_t OldMu100() const {return OldMu100_;}
  Bool_t TkMu100() const {return TkMu100_;}
  Bool_t DiPFJetAve40() const {return DiPFJetAve40_;}
  Bool_t DiPFJetAve60() const {return DiPFJetAve60_;}
  Bool_t DiPFJetAve80() const {return DiPFJetAve80_;}
  Bool_t DiPFJetAve140() const {return DiPFJetAve140_;}
  Bool_t DiPFJetAve200() const {return DiPFJetAve200_;}
  Bool_t DiPFJetAve260() const {return DiPFJetAve260_;}
  Bool_t DiPFJetAve320() const {return DiPFJetAve320_;}
  Bool_t DiPFJetAve400() const {return DiPFJetAve400_;}
  Bool_t DiPFJetAve500() const {return DiPFJetAve500_;}
  Bool_t DiPFJetAve60_HFJEC() const {return DiPFJetAve60_HFJEC_;}
  Bool_t DiPFJetAve80_HFJEC() const {return DiPFJetAve80_HFJEC_;}
  Bool_t DiPFJetAve100_HFJEC() const {return DiPFJetAve100_HFJEC_;}
  Bool_t DiPFJetAve160_HFJEC() const {return DiPFJetAve160_HFJEC_;}
  Bool_t DiPFJetAve220_HFJEC() const {return DiPFJetAve220_HFJEC_;}
  Bool_t DiPFJetAve300_HFJEC() const {return DiPFJetAve300_HFJEC_;}
  Bool_t AK8PFJet15() const {return AK8PFJet15_;}
  Bool_t AK8PFJet25() const {return AK8PFJet25_;}
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
  Bool_t AK8PFJet550() const {return AK8PFJet550_;}
  Bool_t PFJet15() const {return PFJet15_;}
  Bool_t PFJet25() const {return PFJet25_;}
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
  Bool_t PFJet550() const {return PFJet550_;}
  Bool_t PFJetFwd15() const {return PFJetFwd15_;}
  Bool_t PFJetFwd25() const {return PFJetFwd25_;}
  Bool_t PFJetFwd40() const {return PFJetFwd40_;}
  Bool_t PFJetFwd60() const {return PFJetFwd60_;}
  Bool_t PFJetFwd80() const {return PFJetFwd80_;}
  Bool_t PFJetFwd140() const {return PFJetFwd140_;}
  Bool_t PFJetFwd200() const {return PFJetFwd200_;}
  Bool_t PFJetFwd260() const {return PFJetFwd260_;}
  Bool_t PFJetFwd320() const {return PFJetFwd320_;}
  Bool_t PFJetFwd400() const {return PFJetFwd400_;}
  Bool_t PFJetFwd450() const {return PFJetFwd450_;}
  Bool_t PFJetFwd500() const {return PFJetFwd500_;}
  Bool_t AK8PFJetFwd15() const {return AK8PFJetFwd15_;}
  Bool_t AK8PFJetFwd25() const {return AK8PFJetFwd25_;}
  Bool_t AK8PFJetFwd40() const {return AK8PFJetFwd40_;}
  Bool_t AK8PFJetFwd60() const {return AK8PFJetFwd60_;}
  Bool_t AK8PFJetFwd80() const {return AK8PFJetFwd80_;}
  Bool_t AK8PFJetFwd140() const {return AK8PFJetFwd140_;}
  Bool_t AK8PFJetFwd200() const {return AK8PFJetFwd200_;}
  Bool_t AK8PFJetFwd260() const {return AK8PFJetFwd260_;}
  Bool_t AK8PFJetFwd320() const {return AK8PFJetFwd320_;}
  Bool_t AK8PFJetFwd400() const {return AK8PFJetFwd400_;}
  Bool_t AK8PFJetFwd450() const {return AK8PFJetFwd450_;}
  Bool_t AK8PFJetFwd500() const {return AK8PFJetFwd500_;}
  Bool_t PFHT180() const {return PFHT180_;}
  Bool_t PFHT250() const {return PFHT250_;}
  Bool_t PFHT370() const {return PFHT370_;}
  Bool_t PFHT430() const {return PFHT430_;}
  Bool_t PFHT510() const {return PFHT510_;}
  Bool_t PFHT590() const {return PFHT590_;}
  Bool_t PFHT680() const {return PFHT680_;}
  Bool_t PFHT780() const {return PFHT780_;}
  Bool_t PFHT890() const {return PFHT890_;}
  Bool_t PFHT1050() const {return PFHT1050_;}
  Bool_t PFHT500_PFMET100_PFMHT100_IDTight() const {return PFHT500_PFMET100_PFMHT100_IDTight_;}
  Bool_t PFHT500_PFMET110_PFMHT110_IDTight() const {return PFHT500_PFMET110_PFMHT110_IDTight_;}
  Bool_t PFHT700_PFMET85_PFMHT85_IDTight() const {return PFHT700_PFMET85_PFMHT85_IDTight_;}
  Bool_t PFHT700_PFMET95_PFMHT95_IDTight() const {return PFHT700_PFMET95_PFMHT95_IDTight_;}
  Bool_t PFHT800_PFMET75_PFMHT75_IDTight() const {return PFHT800_PFMET75_PFMHT75_IDTight_;}
  Bool_t PFHT800_PFMET85_PFMHT85_IDTight() const {return PFHT800_PFMET85_PFMHT85_IDTight_;}
  Bool_t PFMET110_PFMHT110_IDTight() const {return PFMET110_PFMHT110_IDTight_;}
  Bool_t PFMET120_PFMHT120_IDTight() const {return PFMET120_PFMHT120_IDTight_;}
  Bool_t PFMET130_PFMHT130_IDTight() const {return PFMET130_PFMHT130_IDTight_;}
  Bool_t PFMET140_PFMHT140_IDTight() const {return PFMET140_PFMHT140_IDTight_;}
  Bool_t PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1() const {return PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_;}
  Bool_t PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1() const {return PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_;}
  Bool_t PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1() const {return PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_;}
  Bool_t PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1() const {return PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_;}
  Bool_t PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1() const {return PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_;}
  Bool_t PFMET120_PFMHT120_IDTight_PFHT60() const {return PFMET120_PFMHT120_IDTight_PFHT60_;}
  Bool_t PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60() const {return PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_;}
  Bool_t PFMETTypeOne120_PFMHT120_IDTight_PFHT60() const {return PFMETTypeOne120_PFMHT120_IDTight_PFHT60_;}
  Bool_t PFMETTypeOne110_PFMHT110_IDTight() const {return PFMETTypeOne110_PFMHT110_IDTight_;}
  Bool_t PFMETTypeOne120_PFMHT120_IDTight() const {return PFMETTypeOne120_PFMHT120_IDTight_;}
  Bool_t PFMETTypeOne130_PFMHT130_IDTight() const {return PFMETTypeOne130_PFMHT130_IDTight_;}
  Bool_t PFMETTypeOne140_PFMHT140_IDTight() const {return PFMETTypeOne140_PFMHT140_IDTight_;}
  Bool_t PFMETNoMu110_PFMHTNoMu110_IDTight() const {return PFMETNoMu110_PFMHTNoMu110_IDTight_;}
  Bool_t PFMETNoMu120_PFMHTNoMu120_IDTight() const {return PFMETNoMu120_PFMHTNoMu120_IDTight_;}
  Bool_t PFMETNoMu130_PFMHTNoMu130_IDTight() const {return PFMETNoMu130_PFMHTNoMu130_IDTight_;}
  Bool_t PFMETNoMu140_PFMHTNoMu140_IDTight() const {return PFMETNoMu140_PFMHTNoMu140_IDTight_;}
  Bool_t MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight() const {return MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_;}
  Bool_t MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight() const {return MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_;}
  Bool_t MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight() const {return MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_;}
  Bool_t MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight() const {return MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_;}
  Bool_t L1ETMHadSeeds() const {return L1ETMHadSeeds_;}
  Bool_t CaloMHT90() const {return CaloMHT90_;}
  Bool_t CaloMET80_NotCleaned() const {return CaloMET80_NotCleaned_;}
  Bool_t CaloMET90_NotCleaned() const {return CaloMET90_NotCleaned_;}
  Bool_t CaloMET100_NotCleaned() const {return CaloMET100_NotCleaned_;}
  Bool_t CaloMET110_NotCleaned() const {return CaloMET110_NotCleaned_;}
  Bool_t CaloMET250_NotCleaned() const {return CaloMET250_NotCleaned_;}
  Bool_t CaloMET70_HBHECleaned() const {return CaloMET70_HBHECleaned_;}
  Bool_t CaloMET80_HBHECleaned() const {return CaloMET80_HBHECleaned_;}
  Bool_t CaloMET90_HBHECleaned() const {return CaloMET90_HBHECleaned_;}
  Bool_t CaloMET100_HBHECleaned() const {return CaloMET100_HBHECleaned_;}
  Bool_t CaloMET250_HBHECleaned() const {return CaloMET250_HBHECleaned_;}
  Bool_t CaloMET300_HBHECleaned() const {return CaloMET300_HBHECleaned_;}
  Bool_t CaloMET350_HBHECleaned() const {return CaloMET350_HBHECleaned_;}
  Bool_t PFMET200_NotCleaned() const {return PFMET200_NotCleaned_;}
  Bool_t PFMET200_HBHECleaned() const {return PFMET200_HBHECleaned_;}
  Bool_t PFMET250_HBHECleaned() const {return PFMET250_HBHECleaned_;}
  Bool_t PFMET300_HBHECleaned() const {return PFMET300_HBHECleaned_;}
  Bool_t PFMET200_HBHE_BeamHaloCleaned() const {return PFMET200_HBHE_BeamHaloCleaned_;}
  Bool_t PFMETTypeOne200_HBHE_BeamHaloCleaned() const {return PFMETTypeOne200_HBHE_BeamHaloCleaned_;}
  Bool_t MET105_IsoTrk50() const {return MET105_IsoTrk50_;}
  Bool_t MET120_IsoTrk50() const {return MET120_IsoTrk50_;}
  Bool_t SingleJet30_Mu12_SinglePFJet40() const {return SingleJet30_Mu12_SinglePFJet40_;}
  Bool_t Mu12_DoublePFJets40_CaloBTagDeepCSV_p71() const {return Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_;}
  Bool_t Mu12_DoublePFJets100_CaloBTagDeepCSV_p71() const {return Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_;}
  Bool_t Mu12_DoublePFJets200_CaloBTagDeepCSV_p71() const {return Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_;}
  Bool_t Mu12_DoublePFJets350_CaloBTagDeepCSV_p71() const {return Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_;}
  Bool_t Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71() const {return Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;}
  Bool_t Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71() const {return Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;}
  Bool_t Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71() const {return Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;}
  Bool_t DoublePFJets40_CaloBTagDeepCSV_p71() const {return DoublePFJets40_CaloBTagDeepCSV_p71_;}
  Bool_t DoublePFJets100_CaloBTagDeepCSV_p71() const {return DoublePFJets100_CaloBTagDeepCSV_p71_;}
  Bool_t DoublePFJets200_CaloBTagDeepCSV_p71() const {return DoublePFJets200_CaloBTagDeepCSV_p71_;}
  Bool_t DoublePFJets350_CaloBTagDeepCSV_p71() const {return DoublePFJets350_CaloBTagDeepCSV_p71_;}
  Bool_t DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71() const {return DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;}
  Bool_t DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71() const {return DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;}
  Bool_t Photon300_NoHE() const {return Photon300_NoHE_;}
  Bool_t Mu8_TrkIsoVVL() const {return Mu8_TrkIsoVVL_;}
  Bool_t Mu8_DiEle12_CaloIdL_TrackIdL_DZ() const {return Mu8_DiEle12_CaloIdL_TrackIdL_DZ_;}
  Bool_t Mu8_DiEle12_CaloIdL_TrackIdL() const {return Mu8_DiEle12_CaloIdL_TrackIdL_;}
  Bool_t Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ() const {return Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_;}
  Bool_t Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350() const {return Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_;}
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ() const {return Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_;}
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30() const {return Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_;}
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30() const {return Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_;}
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5() const {return Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_;}
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5() const {return Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_;}
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL() const {return Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_;}
  Bool_t Mu17_TrkIsoVVL() const {return Mu17_TrkIsoVVL_;}
  Bool_t Mu19_TrkIsoVVL() const {return Mu19_TrkIsoVVL_;}
  Bool_t BTagMu_AK4DiJet20_Mu5() const {return BTagMu_AK4DiJet20_Mu5_;}
  Bool_t BTagMu_AK4DiJet40_Mu5() const {return BTagMu_AK4DiJet40_Mu5_;}
  Bool_t BTagMu_AK4DiJet70_Mu5() const {return BTagMu_AK4DiJet70_Mu5_;}
  Bool_t BTagMu_AK4DiJet110_Mu5() const {return BTagMu_AK4DiJet110_Mu5_;}
  Bool_t BTagMu_AK4DiJet170_Mu5() const {return BTagMu_AK4DiJet170_Mu5_;}
  Bool_t BTagMu_AK4Jet300_Mu5() const {return BTagMu_AK4Jet300_Mu5_;}
  Bool_t BTagMu_AK8DiJet170_Mu5() const {return BTagMu_AK8DiJet170_Mu5_;}
  Bool_t BTagMu_AK8Jet170_DoubleMu5() const {return BTagMu_AK8Jet170_DoubleMu5_;}
  Bool_t BTagMu_AK8Jet300_Mu5() const {return BTagMu_AK8Jet300_Mu5_;}
  Bool_t BTagMu_AK4DiJet20_Mu5_noalgo() const {return BTagMu_AK4DiJet20_Mu5_noalgo_;}
  Bool_t BTagMu_AK4DiJet40_Mu5_noalgo() const {return BTagMu_AK4DiJet40_Mu5_noalgo_;}
  Bool_t BTagMu_AK4DiJet70_Mu5_noalgo() const {return BTagMu_AK4DiJet70_Mu5_noalgo_;}
  Bool_t BTagMu_AK4DiJet110_Mu5_noalgo() const {return BTagMu_AK4DiJet110_Mu5_noalgo_;}
  Bool_t BTagMu_AK4DiJet170_Mu5_noalgo() const {return BTagMu_AK4DiJet170_Mu5_noalgo_;}
  Bool_t BTagMu_AK4Jet300_Mu5_noalgo() const {return BTagMu_AK4Jet300_Mu5_noalgo_;}
  Bool_t BTagMu_AK8DiJet170_Mu5_noalgo() const {return BTagMu_AK8DiJet170_Mu5_noalgo_;}
  Bool_t BTagMu_AK8Jet170_DoubleMu5_noalgo() const {return BTagMu_AK8Jet170_DoubleMu5_noalgo_;}
  Bool_t BTagMu_AK8Jet300_Mu5_noalgo() const {return BTagMu_AK8Jet300_Mu5_noalgo_;}
  Bool_t Ele15_Ele8_CaloIdL_TrackIdL_IsoVL() const {return Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_;}
  Bool_t Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() const {return Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;}
  Bool_t Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() const {return Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_;}
  Bool_t Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() const {return Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;}
  Bool_t Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL() const {return Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_;}
  Bool_t Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL() const {return Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_;}
  Bool_t Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ() const {return Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_;}
  Bool_t Mu12_DoublePhoton20() const {return Mu12_DoublePhoton20_;}
  Bool_t TriplePhoton_20_20_20_CaloIdLV2() const {return TriplePhoton_20_20_20_CaloIdLV2_;}
  Bool_t TriplePhoton_20_20_20_CaloIdLV2_R9IdVL() const {return TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_;}
  Bool_t TriplePhoton_30_30_10_CaloIdLV2() const {return TriplePhoton_30_30_10_CaloIdLV2_;}
  Bool_t TriplePhoton_30_30_10_CaloIdLV2_R9IdVL() const {return TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_;}
  Bool_t TriplePhoton_35_35_5_CaloIdLV2_R9IdVL() const {return TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_;}
  Bool_t Photon20() const {return Photon20_;}
  Bool_t Photon33() const {return Photon33_;}
  Bool_t Photon50() const {return Photon50_;}
  Bool_t Photon75() const {return Photon75_;}
  Bool_t Photon90() const {return Photon90_;}
  Bool_t Photon120() const {return Photon120_;}
  Bool_t Photon150() const {return Photon150_;}
  Bool_t Photon175() const {return Photon175_;}
  Bool_t Photon200() const {return Photon200_;}
  Bool_t Photon100EB_TightID_TightIso() const {return Photon100EB_TightID_TightIso_;}
  Bool_t Photon110EB_TightID_TightIso() const {return Photon110EB_TightID_TightIso_;}
  Bool_t Photon120EB_TightID_TightIso() const {return Photon120EB_TightID_TightIso_;}
  Bool_t Photon100EBHE10() const {return Photon100EBHE10_;}
  Bool_t Photon100EEHE10() const {return Photon100EEHE10_;}
  Bool_t Photon100EE_TightID_TightIso() const {return Photon100EE_TightID_TightIso_;}
  Bool_t Photon50_R9Id90_HE10_IsoM() const {return Photon50_R9Id90_HE10_IsoM_;}
  Bool_t Photon75_R9Id90_HE10_IsoM() const {return Photon75_R9Id90_HE10_IsoM_;}
  Bool_t Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3() const {return Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_;}
  Bool_t Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3() const {return Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_;}
  Bool_t Photon90_R9Id90_HE10_IsoM() const {return Photon90_R9Id90_HE10_IsoM_;}
  Bool_t Photon120_R9Id90_HE10_IsoM() const {return Photon120_R9Id90_HE10_IsoM_;}
  Bool_t Photon165_R9Id90_HE10_IsoM() const {return Photon165_R9Id90_HE10_IsoM_;}
  Bool_t Photon90_CaloIdL_PFHT700() const {return Photon90_CaloIdL_PFHT700_;}
  Bool_t Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() const {return Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_;}
  Bool_t Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95() const {return Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_;}
  Bool_t Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55() const {return Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_;}
  Bool_t Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55() const {return Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_;}
  Bool_t Photon35_TwoProngs35() const {return Photon35_TwoProngs35_;}
  Bool_t IsoMu24_TwoProngs35() const {return IsoMu24_TwoProngs35_;}
  Bool_t Dimuon0_Jpsi_L1_NoOS() const {return Dimuon0_Jpsi_L1_NoOS_;}
  Bool_t Dimuon0_Jpsi_NoVertexing_NoOS() const {return Dimuon0_Jpsi_NoVertexing_NoOS_;}
  Bool_t Dimuon0_Jpsi() const {return Dimuon0_Jpsi_;}
  Bool_t Dimuon0_Jpsi_NoVertexing() const {return Dimuon0_Jpsi_NoVertexing_;}
  Bool_t Dimuon0_Jpsi_L1_4R_0er1p5R() const {return Dimuon0_Jpsi_L1_4R_0er1p5R_;}
  Bool_t Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R() const {return Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_;}
  Bool_t Dimuon0_Jpsi3p5_Muon2() const {return Dimuon0_Jpsi3p5_Muon2_;}
  Bool_t Dimuon0_Upsilon_L1_4p5() const {return Dimuon0_Upsilon_L1_4p5_;}
  Bool_t Dimuon0_Upsilon_L1_5() const {return Dimuon0_Upsilon_L1_5_;}
  Bool_t Dimuon0_Upsilon_L1_4p5NoOS() const {return Dimuon0_Upsilon_L1_4p5NoOS_;}
  Bool_t Dimuon0_Upsilon_L1_4p5er2p0() const {return Dimuon0_Upsilon_L1_4p5er2p0_;}
  Bool_t Dimuon0_Upsilon_L1_4p5er2p0M() const {return Dimuon0_Upsilon_L1_4p5er2p0M_;}
  Bool_t Dimuon0_Upsilon_NoVertexing() const {return Dimuon0_Upsilon_NoVertexing_;}
  Bool_t Dimuon0_Upsilon_L1_5M() const {return Dimuon0_Upsilon_L1_5M_;}
  Bool_t Dimuon0_LowMass_L1_0er1p5R() const {return Dimuon0_LowMass_L1_0er1p5R_;}
  Bool_t Dimuon0_LowMass_L1_0er1p5() const {return Dimuon0_LowMass_L1_0er1p5_;}
  Bool_t Dimuon0_LowMass() const {return Dimuon0_LowMass_;}
  Bool_t Dimuon0_LowMass_L1_4() const {return Dimuon0_LowMass_L1_4_;}
  Bool_t Dimuon0_LowMass_L1_4R() const {return Dimuon0_LowMass_L1_4R_;}
  Bool_t Dimuon0_LowMass_L1_TM530() const {return Dimuon0_LowMass_L1_TM530_;}
  Bool_t Dimuon0_Upsilon_Muon_L1_TM0() const {return Dimuon0_Upsilon_Muon_L1_TM0_;}
  Bool_t Dimuon0_Upsilon_Muon_NoL1Mass() const {return Dimuon0_Upsilon_Muon_NoL1Mass_;}
  Bool_t TripleMu_5_3_3_Mass3p8_DZ() const {return TripleMu_5_3_3_Mass3p8_DZ_;}
  Bool_t TripleMu_10_5_5_DZ() const {return TripleMu_10_5_5_DZ_;}
  Bool_t TripleMu_12_10_5() const {return TripleMu_12_10_5_;}
  Bool_t Tau3Mu_Mu7_Mu1_TkMu1_Tau15() const {return Tau3Mu_Mu7_Mu1_TkMu1_Tau15_;}
  Bool_t Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1() const {return Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_;}
  Bool_t Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15() const {return Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_;}
  Bool_t Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1() const {return Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_;}
  Bool_t DoubleMu3_DZ_PFMET50_PFMHT60() const {return DoubleMu3_DZ_PFMET50_PFMHT60_;}
  Bool_t DoubleMu3_DZ_PFMET70_PFMHT70() const {return DoubleMu3_DZ_PFMET70_PFMHT70_;}
  Bool_t DoubleMu3_DZ_PFMET90_PFMHT90() const {return DoubleMu3_DZ_PFMET90_PFMHT90_;}
  Bool_t DoubleMu3_Trk_Tau3mu_NoL1Mass() const {return DoubleMu3_Trk_Tau3mu_NoL1Mass_;}
  Bool_t DoubleMu4_Jpsi_Displaced() const {return DoubleMu4_Jpsi_Displaced_;}
  Bool_t DoubleMu4_Jpsi_NoVertexing() const {return DoubleMu4_Jpsi_NoVertexing_;}
  Bool_t DoubleMu4_JpsiTrkTrk_Displaced() const {return DoubleMu4_JpsiTrkTrk_Displaced_;}
  Bool_t DoubleMu43NoFiltersNoVtx() const {return DoubleMu43NoFiltersNoVtx_;}
  Bool_t DoubleMu48NoFiltersNoVtx() const {return DoubleMu48NoFiltersNoVtx_;}
  Bool_t Mu43NoFiltersNoVtx_Photon43_CaloIdL() const {return Mu43NoFiltersNoVtx_Photon43_CaloIdL_;}
  Bool_t Mu48NoFiltersNoVtx_Photon48_CaloIdL() const {return Mu48NoFiltersNoVtx_Photon48_CaloIdL_;}
  Bool_t Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL() const {return Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_;}
  Bool_t Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL() const {return Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_;}
  Bool_t DoubleMu33NoFiltersNoVtxDisplaced() const {return DoubleMu33NoFiltersNoVtxDisplaced_;}
  Bool_t DoubleMu40NoFiltersNoVtxDisplaced() const {return DoubleMu40NoFiltersNoVtxDisplaced_;}
  Bool_t DoubleMu20_7_Mass0to30_L1_DM4() const {return DoubleMu20_7_Mass0to30_L1_DM4_;}
  Bool_t DoubleMu20_7_Mass0to30_L1_DM4EG() const {return DoubleMu20_7_Mass0to30_L1_DM4EG_;}
  Bool_t HT425() const {return HT425_;}
  Bool_t HT430_DisplacedDijet40_DisplacedTrack() const {return HT430_DisplacedDijet40_DisplacedTrack_;}
  Bool_t HT500_DisplacedDijet40_DisplacedTrack() const {return HT500_DisplacedDijet40_DisplacedTrack_;}
  Bool_t HT430_DisplacedDijet60_DisplacedTrack() const {return HT430_DisplacedDijet60_DisplacedTrack_;}
  Bool_t HT400_DisplacedDijet40_DisplacedTrack() const {return HT400_DisplacedDijet40_DisplacedTrack_;}
  Bool_t HT650_DisplacedDijet60_Inclusive() const {return HT650_DisplacedDijet60_Inclusive_;}
  Bool_t HT550_DisplacedDijet60_Inclusive() const {return HT550_DisplacedDijet60_Inclusive_;}
  Bool_t DiJet110_35_Mjj650_PFMET110() const {return DiJet110_35_Mjj650_PFMET110_;}
  Bool_t DiJet110_35_Mjj650_PFMET120() const {return DiJet110_35_Mjj650_PFMET120_;}
  Bool_t DiJet110_35_Mjj650_PFMET130() const {return DiJet110_35_Mjj650_PFMET130_;}
  Bool_t TripleJet110_35_35_Mjj650_PFMET110() const {return TripleJet110_35_35_Mjj650_PFMET110_;}
  Bool_t TripleJet110_35_35_Mjj650_PFMET120() const {return TripleJet110_35_35_Mjj650_PFMET120_;}
  Bool_t TripleJet110_35_35_Mjj650_PFMET130() const {return TripleJet110_35_35_Mjj650_PFMET130_;}
  Bool_t Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned() const {return Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_;}
  Bool_t Ele28_eta2p1_WPTight_Gsf_HT150() const {return Ele28_eta2p1_WPTight_Gsf_HT150_;}
  Bool_t Ele28_HighEta_SC20_Mass55() const {return Ele28_HighEta_SC20_Mass55_;}
  Bool_t DoubleMu20_7_Mass0to30_Photon23() const {return DoubleMu20_7_Mass0to30_Photon23_;}
  Bool_t Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5() const {return Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_;}
  Bool_t Ele15_IsoVVVL_PFHT450_PFMET50() const {return Ele15_IsoVVVL_PFHT450_PFMET50_;}
  Bool_t Ele15_IsoVVVL_PFHT450() const {return Ele15_IsoVVVL_PFHT450_;}
  Bool_t Ele50_IsoVVVL_PFHT450() const {return Ele50_IsoVVVL_PFHT450_;}
  Bool_t Ele15_IsoVVVL_PFHT600() const {return Ele15_IsoVVVL_PFHT600_;}
  Bool_t Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60() const {return Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_;}
  Bool_t Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60() const {return Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_;}
  Bool_t Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60() const {return Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_;}
  Bool_t Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5() const {return Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_;}
  Bool_t Mu15_IsoVVVL_PFHT450_PFMET50() const {return Mu15_IsoVVVL_PFHT450_PFMET50_;}
  Bool_t Mu15_IsoVVVL_PFHT450() const {return Mu15_IsoVVVL_PFHT450_;}
  Bool_t Mu50_IsoVVVL_PFHT450() const {return Mu50_IsoVVVL_PFHT450_;}
  Bool_t Mu15_IsoVVVL_PFHT600() const {return Mu15_IsoVVVL_PFHT600_;}
  Bool_t Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight() const {return Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_;}
  Bool_t Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight() const {return Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_;}
  Bool_t Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight() const {return Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_;}
  Bool_t Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight() const {return Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_;}
  Bool_t Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight() const {return Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_;}
  Bool_t Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight() const {return Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_;}
  Bool_t Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight() const {return Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_;}
  Bool_t Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight() const {return Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_;}
  Bool_t Dimuon10_PsiPrime_Barrel_Seagulls() const {return Dimuon10_PsiPrime_Barrel_Seagulls_;}
  Bool_t Dimuon20_Jpsi_Barrel_Seagulls() const {return Dimuon20_Jpsi_Barrel_Seagulls_;}
  Bool_t Dimuon12_Upsilon_y1p4() const {return Dimuon12_Upsilon_y1p4_;}
  Bool_t Dimuon14_Phi_Barrel_Seagulls() const {return Dimuon14_Phi_Barrel_Seagulls_;}
  Bool_t Dimuon18_PsiPrime() const {return Dimuon18_PsiPrime_;}
  Bool_t Dimuon25_Jpsi() const {return Dimuon25_Jpsi_;}
  Bool_t Dimuon18_PsiPrime_noCorrL1() const {return Dimuon18_PsiPrime_noCorrL1_;}
  Bool_t Dimuon24_Upsilon_noCorrL1() const {return Dimuon24_Upsilon_noCorrL1_;}
  Bool_t Dimuon24_Phi_noCorrL1() const {return Dimuon24_Phi_noCorrL1_;}
  Bool_t Dimuon25_Jpsi_noCorrL1() const {return Dimuon25_Jpsi_noCorrL1_;}
  Bool_t DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8() const {return DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_;}
  Bool_t DiMu9_Ele9_CaloIdL_TrackIdL_DZ() const {return DiMu9_Ele9_CaloIdL_TrackIdL_DZ_;}
  Bool_t DiMu9_Ele9_CaloIdL_TrackIdL() const {return DiMu9_Ele9_CaloIdL_TrackIdL_;}
  Bool_t DoubleIsoMu20_eta2p1() const {return DoubleIsoMu20_eta2p1_;}
  Bool_t TrkMu12_DoubleTrkMu5NoFiltersNoVtx() const {return TrkMu12_DoubleTrkMu5NoFiltersNoVtx_;}
  Bool_t TrkMu16_DoubleTrkMu6NoFiltersNoVtx() const {return TrkMu16_DoubleTrkMu6NoFiltersNoVtx_;}
  Bool_t TrkMu17_DoubleTrkMu8NoFiltersNoVtx() const {return TrkMu17_DoubleTrkMu8NoFiltersNoVtx_;}
  Bool_t Mu8() const {return Mu8_;}
  Bool_t Mu17() const {return Mu17_;}
  Bool_t Mu19() const {return Mu19_;}
  Bool_t Mu17_Photon30_IsoCaloId() const {return Mu17_Photon30_IsoCaloId_;}
  Bool_t Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30() const {return Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_;}
  Bool_t Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30() const {return Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_;}
  Bool_t Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30() const {return Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_;}
  Bool_t Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30() const {return Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_;}
  Bool_t Ele8_CaloIdM_TrackIdM_PFJet30() const {return Ele8_CaloIdM_TrackIdM_PFJet30_;}
  Bool_t Ele17_CaloIdM_TrackIdM_PFJet30() const {return Ele17_CaloIdM_TrackIdM_PFJet30_;}
  Bool_t Ele23_CaloIdM_TrackIdM_PFJet30() const {return Ele23_CaloIdM_TrackIdM_PFJet30_;}
  Bool_t Ele50_CaloIdVT_GsfTrkIdT_PFJet165() const {return Ele50_CaloIdVT_GsfTrkIdT_PFJet165_;}
  Bool_t Ele115_CaloIdVT_GsfTrkIdT() const {return Ele115_CaloIdVT_GsfTrkIdT_;}
  Bool_t Ele135_CaloIdVT_GsfTrkIdT() const {return Ele135_CaloIdVT_GsfTrkIdT_;}
  Bool_t Ele145_CaloIdVT_GsfTrkIdT() const {return Ele145_CaloIdVT_GsfTrkIdT_;}
  Bool_t Ele200_CaloIdVT_GsfTrkIdT() const {return Ele200_CaloIdVT_GsfTrkIdT_;}
  Bool_t Ele250_CaloIdVT_GsfTrkIdT() const {return Ele250_CaloIdVT_GsfTrkIdT_;}
  Bool_t Ele300_CaloIdVT_GsfTrkIdT() const {return Ele300_CaloIdVT_GsfTrkIdT_;}
  Bool_t PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5() const {return PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_;}
  Bool_t PFHT330PT30_QuadPFJet_75_60_45_40() const {return PFHT330PT30_QuadPFJet_75_60_45_40_;}
  Bool_t PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94() const {return PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_;}
  Bool_t PFHT400_SixPFJet32() const {return PFHT400_SixPFJet32_;}
  Bool_t PFHT450_SixPFJet36_PFBTagDeepCSV_1p59() const {return PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_;}
  Bool_t PFHT450_SixPFJet36() const {return PFHT450_SixPFJet36_;}
  Bool_t PFHT350() const {return PFHT350_;}
  Bool_t PFHT350MinPFJet15() const {return PFHT350MinPFJet15_;}
  Bool_t Photon60_R9Id90_CaloIdL_IsoL() const {return Photon60_R9Id90_CaloIdL_IsoL_;}
  Bool_t Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL() const {return Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_;}
  Bool_t Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15() const {return Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_;}
  Bool_t ECALHT800() const {return ECALHT800_;}
  Bool_t DiSC30_18_EIso_AND_HE_Mass70() const {return DiSC30_18_EIso_AND_HE_Mass70_;}
  Bool_t Physics() const {return Physics_;}
  Bool_t Physics_part0() const {return Physics_part0_;}
  Bool_t Physics_part1() const {return Physics_part1_;}
  Bool_t Physics_part2() const {return Physics_part2_;}
  Bool_t Physics_part3() const {return Physics_part3_;}
  Bool_t Physics_part4() const {return Physics_part4_;}
  Bool_t Physics_part5() const {return Physics_part5_;}
  Bool_t Physics_part6() const {return Physics_part6_;}
  Bool_t Physics_part7() const {return Physics_part7_;}
  Bool_t Random() const {return Random_;}
  Bool_t ZeroBias() const {return ZeroBias_;}
  Bool_t ZeroBias_Alignment() const {return ZeroBias_Alignment_;}
  Bool_t ZeroBias_part0() const {return ZeroBias_part0_;}
  Bool_t ZeroBias_part1() const {return ZeroBias_part1_;}
  Bool_t ZeroBias_part2() const {return ZeroBias_part2_;}
  Bool_t ZeroBias_part3() const {return ZeroBias_part3_;}
  Bool_t ZeroBias_part4() const {return ZeroBias_part4_;}
  Bool_t ZeroBias_part5() const {return ZeroBias_part5_;}
  Bool_t ZeroBias_part6() const {return ZeroBias_part6_;}
  Bool_t ZeroBias_part7() const {return ZeroBias_part7_;}
  Bool_t AK4CaloJet30() const {return AK4CaloJet30_;}
  Bool_t AK4CaloJet40() const {return AK4CaloJet40_;}
  Bool_t AK4CaloJet50() const {return AK4CaloJet50_;}
  Bool_t AK4CaloJet80() const {return AK4CaloJet80_;}
  Bool_t AK4CaloJet100() const {return AK4CaloJet100_;}
  Bool_t AK4CaloJet120() const {return AK4CaloJet120_;}
  Bool_t AK4PFJet30() const {return AK4PFJet30_;}
  Bool_t AK4PFJet50() const {return AK4PFJet50_;}
  Bool_t AK4PFJet80() const {return AK4PFJet80_;}
  Bool_t AK4PFJet100() const {return AK4PFJet100_;}
  Bool_t AK4PFJet120() const {return AK4PFJet120_;}
  Bool_t SinglePhoton10_Eta3p1ForPPRef() const {return SinglePhoton10_Eta3p1ForPPRef_;}
  Bool_t SinglePhoton20_Eta3p1ForPPRef() const {return SinglePhoton20_Eta3p1ForPPRef_;}
  Bool_t SinglePhoton30_Eta3p1ForPPRef() const {return SinglePhoton30_Eta3p1ForPPRef_;}
  Bool_t Photon20_HoverELoose() const {return Photon20_HoverELoose_;}
  Bool_t Photon30_HoverELoose() const {return Photon30_HoverELoose_;}
  Bool_t EcalCalibration() const {return EcalCalibration_;}
  Bool_t HcalCalibration() const {return HcalCalibration_;}
  Bool_t L1UnpairedBunchBptxMinus() const {return L1UnpairedBunchBptxMinus_;}
  Bool_t L1UnpairedBunchBptxPlus() const {return L1UnpairedBunchBptxPlus_;}
  Bool_t L1NotBptxOR() const {return L1NotBptxOR_;}
  Bool_t L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142() const {return L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_;}
  Bool_t CDC_L2cosmic_5_er1p0() const {return CDC_L2cosmic_5_er1p0_;}
  Bool_t CDC_L2cosmic_5p5_er1p0() const {return CDC_L2cosmic_5p5_er1p0_;}
  Bool_t HcalNZS() const {return HcalNZS_;}
  Bool_t HcalPhiSym() const {return HcalPhiSym_;}
  Bool_t HcalIsolatedbunch() const {return HcalIsolatedbunch_;}
  Bool_t IsoTrackHB() const {return IsoTrackHB_;}
  Bool_t IsoTrackHE() const {return IsoTrackHE_;}
  Bool_t ZeroBias_FirstCollisionAfterAbortGap() const {return ZeroBias_FirstCollisionAfterAbortGap_;}
  Bool_t ZeroBias_IsolatedBunches() const {return ZeroBias_IsolatedBunches_;}
  Bool_t ZeroBias_FirstCollisionInTrain() const {return ZeroBias_FirstCollisionInTrain_;}
  Bool_t ZeroBias_LastCollisionInTrain() const {return ZeroBias_LastCollisionInTrain_;}
  Bool_t ZeroBias_FirstBXAfterTrain() const {return ZeroBias_FirstBXAfterTrain_;}
  Bool_t IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr() const {return IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_;}
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90() const {return MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_;}
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100() const {return MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_;}
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110() const {return MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_;}
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120() const {return MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_;}
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130() const {return MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_;}
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140() const {return MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_;}
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr() const {return MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_;}
  Bool_t MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr() const {return MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_;}
  Bool_t MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1() const {return MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_;}
  Bool_t MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1() const {return MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_;}
  Bool_t MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1() const {return MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_;}
  Bool_t Ele16_Ele12_Ele8_CaloIdL_TrackIdL() const {return Ele16_Ele12_Ele8_CaloIdL_TrackIdL_;}
  Bool_t Rsq0p35() const {return Rsq0p35_;}
  Bool_t Rsq0p40() const {return Rsq0p40_;}
  Bool_t RsqMR300_Rsq0p09_MR200() const {return RsqMR300_Rsq0p09_MR200_;}
  Bool_t RsqMR320_Rsq0p09_MR200() const {return RsqMR320_Rsq0p09_MR200_;}
  Bool_t RsqMR300_Rsq0p09_MR200_4jet() const {return RsqMR300_Rsq0p09_MR200_4jet_;}
  Bool_t RsqMR320_Rsq0p09_MR200_4jet() const {return RsqMR320_Rsq0p09_MR200_4jet_;}
  Bool_t IsoMu27_MET90() const {return IsoMu27_MET90_;}
  Bool_t DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg() const {return DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_;}
  Bool_t DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg() const {return DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_;}
  Bool_t DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg() const {return DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_;}
  Bool_t DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg() const {return DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_;}
  Bool_t DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg() const {return DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_;}
  Bool_t DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg() const {return DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_;}
  Bool_t DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg() const {return DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_;}
  Bool_t DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg() const {return DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_;}
  Bool_t VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1() const {return VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_;}
  Bool_t VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1() const {return VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_;}
  Bool_t VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1() const {return VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_;}
  Bool_t Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50() const {return Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_;}
  Bool_t Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3() const {return Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_;}
  Bool_t Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3() const {return Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_;}
  Bool_t PFMET100_PFMHT100_IDTight_PFHT60() const {return PFMET100_PFMHT100_IDTight_PFHT60_;}
  Bool_t PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60() const {return PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_;}
  Bool_t PFMETTypeOne100_PFMHT100_IDTight_PFHT60() const {return PFMETTypeOne100_PFMHT100_IDTight_PFHT60_;}
  Bool_t Mu18_Mu9_SameSign() const {return Mu18_Mu9_SameSign_;}
  Bool_t Mu18_Mu9_SameSign_DZ() const {return Mu18_Mu9_SameSign_DZ_;}
  Bool_t Mu18_Mu9() const {return Mu18_Mu9_;}
  Bool_t Mu18_Mu9_DZ() const {return Mu18_Mu9_DZ_;}
  Bool_t Mu20_Mu10_SameSign() const {return Mu20_Mu10_SameSign_;}
  Bool_t Mu20_Mu10_SameSign_DZ() const {return Mu20_Mu10_SameSign_DZ_;}
  Bool_t Mu20_Mu10() const {return Mu20_Mu10_;}
  Bool_t Mu20_Mu10_DZ() const {return Mu20_Mu10_DZ_;}
  Bool_t Mu23_Mu12_SameSign() const {return Mu23_Mu12_SameSign_;}
  Bool_t Mu23_Mu12_SameSign_DZ() const {return Mu23_Mu12_SameSign_DZ_;}
  Bool_t Mu23_Mu12() const {return Mu23_Mu12_;}
  Bool_t Mu23_Mu12_DZ() const {return Mu23_Mu12_DZ_;}
  Bool_t DoubleMu2_Jpsi_DoubleTrk1_Phi1p05() const {return DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_;}
  Bool_t DoubleMu2_Jpsi_DoubleTkMu0_Phi() const {return DoubleMu2_Jpsi_DoubleTkMu0_Phi_;}
  Bool_t DoubleMu3_DCA_PFMET50_PFMHT60() const {return DoubleMu3_DCA_PFMET50_PFMHT60_;}
  Bool_t TripleMu_5_3_3_Mass3p8_DCA() const {return TripleMu_5_3_3_Mass3p8_DCA_;}
  Bool_t QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1() const {return QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;}
  Bool_t QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1() const {return QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;}
  Bool_t QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1() const {return QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;}
  Bool_t QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2() const {return QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_;}
  Bool_t QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2() const {return QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_;}
  Bool_t QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2() const {return QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_;}
  Bool_t QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2() const {return QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_;}
  Bool_t QuadPFJet98_83_71_15() const {return QuadPFJet98_83_71_15_;}
  Bool_t QuadPFJet103_88_75_15() const {return QuadPFJet103_88_75_15_;}
  Bool_t QuadPFJet105_88_76_15() const {return QuadPFJet105_88_76_15_;}
  Bool_t QuadPFJet111_90_80_15() const {return QuadPFJet111_90_80_15_;}
  Bool_t AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17() const {return AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_;}
  Bool_t AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1() const {return AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_;}
  Bool_t AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02() const {return AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_;}
  Bool_t AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2() const {return AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_;}
  Bool_t AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4() const {return AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_;}
  Bool_t Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55() const {return Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_;}
  Bool_t Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto() const {return Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_;}
  Bool_t Mu12_IP6_part0() const {return Mu12_IP6_part0_;}
  Bool_t Mu12_IP6_part1() const {return Mu12_IP6_part1_;}
  Bool_t Mu12_IP6_part2() const {return Mu12_IP6_part2_;}
  Bool_t Mu12_IP6_part3() const {return Mu12_IP6_part3_;}
  Bool_t Mu12_IP6_part4() const {return Mu12_IP6_part4_;}
  Bool_t Mu9_IP5_part0() const {return Mu9_IP5_part0_;}
  Bool_t Mu9_IP5_part1() const {return Mu9_IP5_part1_;}
  Bool_t Mu9_IP5_part2() const {return Mu9_IP5_part2_;}
  Bool_t Mu9_IP5_part3() const {return Mu9_IP5_part3_;}
  Bool_t Mu9_IP5_part4() const {return Mu9_IP5_part4_;}
  Bool_t Mu7_IP4_part0() const {return Mu7_IP4_part0_;}
  Bool_t Mu7_IP4_part1() const {return Mu7_IP4_part1_;}
  Bool_t Mu7_IP4_part2() const {return Mu7_IP4_part2_;}
  Bool_t Mu7_IP4_part3() const {return Mu7_IP4_part3_;}
  Bool_t Mu7_IP4_part4() const {return Mu7_IP4_part4_;}
  Bool_t Mu9_IP4_part0() const {return Mu9_IP4_part0_;}
  Bool_t Mu9_IP4_part1() const {return Mu9_IP4_part1_;}
  Bool_t Mu9_IP4_part2() const {return Mu9_IP4_part2_;}
  Bool_t Mu9_IP4_part3() const {return Mu9_IP4_part3_;}
  Bool_t Mu9_IP4_part4() const {return Mu9_IP4_part4_;}
  Bool_t Mu8_IP5_part0() const {return Mu8_IP5_part0_;}
  Bool_t Mu8_IP5_part1() const {return Mu8_IP5_part1_;}
  Bool_t Mu8_IP5_part2() const {return Mu8_IP5_part2_;}
  Bool_t Mu8_IP5_part3() const {return Mu8_IP5_part3_;}
  Bool_t Mu8_IP5_part4() const {return Mu8_IP5_part4_;}
  Bool_t Mu8_IP6_part0() const {return Mu8_IP6_part0_;}
  Bool_t Mu8_IP6_part1() const {return Mu8_IP6_part1_;}
  Bool_t Mu8_IP6_part2() const {return Mu8_IP6_part2_;}
  Bool_t Mu8_IP6_part3() const {return Mu8_IP6_part3_;}
  Bool_t Mu8_IP6_part4() const {return Mu8_IP6_part4_;}
  Bool_t Mu9_IP6_part0() const {return Mu9_IP6_part0_;}
  Bool_t Mu9_IP6_part1() const {return Mu9_IP6_part1_;}
  Bool_t Mu9_IP6_part2() const {return Mu9_IP6_part2_;}
  Bool_t Mu9_IP6_part3() const {return Mu9_IP6_part3_;}
  Bool_t Mu9_IP6_part4() const {return Mu9_IP6_part4_;}
  Bool_t Mu8_IP3_part0() const {return Mu8_IP3_part0_;}
  Bool_t Mu8_IP3_part1() const {return Mu8_IP3_part1_;}
  Bool_t Mu8_IP3_part2() const {return Mu8_IP3_part2_;}
  Bool_t Mu8_IP3_part3() const {return Mu8_IP3_part3_;}
  Bool_t Mu8_IP3_part4() const {return Mu8_IP3_part4_;}
  Bool_t QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1() const {return QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;}
  Bool_t TrkMu6NoFiltersNoVtx() const {return TrkMu6NoFiltersNoVtx_;}
  Bool_t TrkMu16NoFiltersNoVtx() const {return TrkMu16NoFiltersNoVtx_;}
  Bool_t DoubleTrkMu_16_6_NoFiltersNoVtx() const {return DoubleTrkMu_16_6_NoFiltersNoVtx_;}
  Bool_t HLTriggerFinalPath() const {return HLTriggerFinalPath_;}
private:
  Bool_t HLTriggerFirstPath_;
  Bool_t AK8PFJet360_TrimMass30_;
  Bool_t AK8PFJet380_TrimMass30_;
  Bool_t AK8PFJet400_TrimMass30_;
  Bool_t AK8PFJet420_TrimMass30_;
  Bool_t AK8PFHT750_TrimMass50_;
  Bool_t AK8PFHT800_TrimMass50_;
  Bool_t AK8PFHT850_TrimMass50_;
  Bool_t AK8PFHT900_TrimMass50_;
  Bool_t CaloJet500_NoJetID_;
  Bool_t CaloJet550_NoJetID_;
  Bool_t DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_;
  Bool_t DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_;
  Bool_t Trimuon5_3p5_2_Upsilon_Muon_;
  Bool_t TrimuonOpen_5_3p5_2_Upsilon_Muon_;
  Bool_t DoubleEle25_CaloIdL_MW_;
  Bool_t DoubleEle27_CaloIdL_MW_;
  Bool_t DoubleEle33_CaloIdL_MW_;
  Bool_t DoubleEle24_eta2p1_WPTight_Gsf_;
  Bool_t DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_;
  Bool_t DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_;
  Bool_t Ele27_Ele37_CaloIdL_MW_;
  Bool_t Mu27_Ele37_CaloIdL_MW_;
  Bool_t Mu37_Ele27_CaloIdL_MW_;
  Bool_t Mu37_TkMu27_;
  Bool_t DoubleMu4_3_Bs_;
  Bool_t DoubleMu4_3_Jpsi_;
  Bool_t DoubleMu4_JpsiTrk_Displaced_;
  Bool_t DoubleMu4_LowMassNonResonantTrk_Displaced_;
  Bool_t DoubleMu3_Trk_Tau3mu_;
  Bool_t DoubleMu3_TkMu_DsTau3Mu_;
  Bool_t DoubleMu4_PsiPrimeTrk_Displaced_;
  Bool_t DoubleMu4_Mass3p8_DZ_PFHT350_;
  Bool_t Mu3_PFJet40_;
  Bool_t Mu7p5_L2Mu2_Jpsi_;
  Bool_t Mu7p5_L2Mu2_Upsilon_;
  Bool_t Mu7p5_Track2_Jpsi_;
  Bool_t Mu7p5_Track3p5_Jpsi_;
  Bool_t Mu7p5_Track7_Jpsi_;
  Bool_t Mu7p5_Track2_Upsilon_;
  Bool_t Mu7p5_Track3p5_Upsilon_;
  Bool_t Mu7p5_Track7_Upsilon_;
  Bool_t Mu3_L1SingleMu5orSingleMu7_;
  Bool_t DoublePhoton33_CaloIdL_;
  Bool_t DoublePhoton70_;
  Bool_t DoublePhoton85_;
  Bool_t Ele20_WPTight_Gsf_;
  Bool_t Ele15_WPLoose_Gsf_;
  Bool_t Ele17_WPLoose_Gsf_;
  Bool_t Ele20_WPLoose_Gsf_;
  Bool_t Ele20_eta2p1_WPLoose_Gsf_;
  Bool_t DiEle27_WPTightCaloOnly_L1DoubleEG_;
  Bool_t Ele27_WPTight_Gsf_;
  Bool_t Ele28_WPTight_Gsf_;
  Bool_t Ele30_WPTight_Gsf_;
  Bool_t Ele32_WPTight_Gsf_;
  Bool_t Ele35_WPTight_Gsf_;
  Bool_t Ele35_WPTight_Gsf_L1EGMT_;
  Bool_t Ele38_WPTight_Gsf_;
  Bool_t Ele40_WPTight_Gsf_;
  Bool_t Ele32_WPTight_Gsf_L1DoubleEG_;
  Bool_t Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_;
  Bool_t Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_;
  Bool_t Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_;
  Bool_t Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;
  Bool_t Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;
  Bool_t Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;
  Bool_t HT450_Beamspot_;
  Bool_t HT300_Beamspot_;
  Bool_t ZeroBias_Beamspot_;
  Bool_t IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_;
  Bool_t IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_;
  Bool_t IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_;
  Bool_t IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;
  Bool_t IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;
  Bool_t IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;
  Bool_t IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_;
  Bool_t IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_;
  Bool_t IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_;
  Bool_t IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_;
  Bool_t IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;
  Bool_t IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;
  Bool_t IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;
  Bool_t IsoMu20_;
  Bool_t IsoMu24_;
  Bool_t IsoMu24_eta2p1_;
  Bool_t IsoMu27_;
  Bool_t IsoMu30_;
  Bool_t UncorrectedJetE30_NoBPTX_;
  Bool_t UncorrectedJetE30_NoBPTX3BX_;
  Bool_t UncorrectedJetE60_NoBPTX3BX_;
  Bool_t UncorrectedJetE70_NoBPTX3BX_;
  Bool_t L1SingleMu18_;
  Bool_t L1SingleMu25_;
  Bool_t L2Mu10_;
  Bool_t L2Mu10_NoVertex_NoBPTX3BX_;
  Bool_t L2Mu10_NoVertex_NoBPTX_;
  Bool_t L2Mu45_NoVertex_3Sta_NoBPTX3BX_;
  Bool_t L2Mu40_NoVertex_3Sta_NoBPTX3BX_;
  Bool_t L2Mu50_;
  Bool_t L2Mu23NoVtx_2Cha_;
  Bool_t L2Mu23NoVtx_2Cha_CosmicSeed_;
  Bool_t DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_;
  Bool_t DoubleL2Mu30NoVtx_2Cha_Eta2p4_;
  Bool_t DoubleL2Mu50_;
  Bool_t DoubleL2Mu23NoVtx_2Cha_CosmicSeed_;
  Bool_t DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_;
  Bool_t DoubleL2Mu25NoVtx_2Cha_CosmicSeed_;
  Bool_t DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_;
  Bool_t DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_;
  Bool_t DoubleL2Mu23NoVtx_2Cha_;
  Bool_t DoubleL2Mu23NoVtx_2Cha_NoL2Matched_;
  Bool_t DoubleL2Mu25NoVtx_2Cha_;
  Bool_t DoubleL2Mu25NoVtx_2Cha_NoL2Matched_;
  Bool_t DoubleL2Mu25NoVtx_2Cha_Eta2p4_;
  Bool_t Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_;
  Bool_t Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_;
  Bool_t Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_;
  Bool_t Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_;
  Bool_t Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_;
  Bool_t Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_;
  Bool_t Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_;
  Bool_t Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_;
  Bool_t Mu25_TkMu0_Onia_;
  Bool_t Mu30_TkMu0_Psi_;
  Bool_t Mu30_TkMu0_Upsilon_;
  Bool_t Mu20_TkMu0_Phi_;
  Bool_t Mu25_TkMu0_Phi_;
  Bool_t Mu12_;
  Bool_t Mu15_;
  Bool_t Mu20_;
  Bool_t Mu27_;
  Bool_t Mu50_;
  Bool_t Mu55_;
  Bool_t OldMu100_;
  Bool_t TkMu100_;
  Bool_t DiPFJetAve40_;
  Bool_t DiPFJetAve60_;
  Bool_t DiPFJetAve80_;
  Bool_t DiPFJetAve140_;
  Bool_t DiPFJetAve200_;
  Bool_t DiPFJetAve260_;
  Bool_t DiPFJetAve320_;
  Bool_t DiPFJetAve400_;
  Bool_t DiPFJetAve500_;
  Bool_t DiPFJetAve60_HFJEC_;
  Bool_t DiPFJetAve80_HFJEC_;
  Bool_t DiPFJetAve100_HFJEC_;
  Bool_t DiPFJetAve160_HFJEC_;
  Bool_t DiPFJetAve220_HFJEC_;
  Bool_t DiPFJetAve300_HFJEC_;
  Bool_t AK8PFJet15_;
  Bool_t AK8PFJet25_;
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
  Bool_t AK8PFJet550_;
  Bool_t PFJet15_;
  Bool_t PFJet25_;
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
  Bool_t PFJet550_;
  Bool_t PFJetFwd15_;
  Bool_t PFJetFwd25_;
  Bool_t PFJetFwd40_;
  Bool_t PFJetFwd60_;
  Bool_t PFJetFwd80_;
  Bool_t PFJetFwd140_;
  Bool_t PFJetFwd200_;
  Bool_t PFJetFwd260_;
  Bool_t PFJetFwd320_;
  Bool_t PFJetFwd400_;
  Bool_t PFJetFwd450_;
  Bool_t PFJetFwd500_;
  Bool_t AK8PFJetFwd15_;
  Bool_t AK8PFJetFwd25_;
  Bool_t AK8PFJetFwd40_;
  Bool_t AK8PFJetFwd60_;
  Bool_t AK8PFJetFwd80_;
  Bool_t AK8PFJetFwd140_;
  Bool_t AK8PFJetFwd200_;
  Bool_t AK8PFJetFwd260_;
  Bool_t AK8PFJetFwd320_;
  Bool_t AK8PFJetFwd400_;
  Bool_t AK8PFJetFwd450_;
  Bool_t AK8PFJetFwd500_;
  Bool_t PFHT180_;
  Bool_t PFHT250_;
  Bool_t PFHT370_;
  Bool_t PFHT430_;
  Bool_t PFHT510_;
  Bool_t PFHT590_;
  Bool_t PFHT680_;
  Bool_t PFHT780_;
  Bool_t PFHT890_;
  Bool_t PFHT1050_;
  Bool_t PFHT500_PFMET100_PFMHT100_IDTight_;
  Bool_t PFHT500_PFMET110_PFMHT110_IDTight_;
  Bool_t PFHT700_PFMET85_PFMHT85_IDTight_;
  Bool_t PFHT700_PFMET95_PFMHT95_IDTight_;
  Bool_t PFHT800_PFMET75_PFMHT75_IDTight_;
  Bool_t PFHT800_PFMET85_PFMHT85_IDTight_;
  Bool_t PFMET110_PFMHT110_IDTight_;
  Bool_t PFMET120_PFMHT120_IDTight_;
  Bool_t PFMET130_PFMHT130_IDTight_;
  Bool_t PFMET140_PFMHT140_IDTight_;
  Bool_t PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t PFMET120_PFMHT120_IDTight_PFHT60_;
  Bool_t PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_;
  Bool_t PFMETTypeOne120_PFMHT120_IDTight_PFHT60_;
  Bool_t PFMETTypeOne110_PFMHT110_IDTight_;
  Bool_t PFMETTypeOne120_PFMHT120_IDTight_;
  Bool_t PFMETTypeOne130_PFMHT130_IDTight_;
  Bool_t PFMETTypeOne140_PFMHT140_IDTight_;
  Bool_t PFMETNoMu110_PFMHTNoMu110_IDTight_;
  Bool_t PFMETNoMu120_PFMHTNoMu120_IDTight_;
  Bool_t PFMETNoMu130_PFMHTNoMu130_IDTight_;
  Bool_t PFMETNoMu140_PFMHTNoMu140_IDTight_;
  Bool_t MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_;
  Bool_t MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_;
  Bool_t MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_;
  Bool_t MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_;
  Bool_t L1ETMHadSeeds_;
  Bool_t CaloMHT90_;
  Bool_t CaloMET80_NotCleaned_;
  Bool_t CaloMET90_NotCleaned_;
  Bool_t CaloMET100_NotCleaned_;
  Bool_t CaloMET110_NotCleaned_;
  Bool_t CaloMET250_NotCleaned_;
  Bool_t CaloMET70_HBHECleaned_;
  Bool_t CaloMET80_HBHECleaned_;
  Bool_t CaloMET90_HBHECleaned_;
  Bool_t CaloMET100_HBHECleaned_;
  Bool_t CaloMET250_HBHECleaned_;
  Bool_t CaloMET300_HBHECleaned_;
  Bool_t CaloMET350_HBHECleaned_;
  Bool_t PFMET200_NotCleaned_;
  Bool_t PFMET200_HBHECleaned_;
  Bool_t PFMET250_HBHECleaned_;
  Bool_t PFMET300_HBHECleaned_;
  Bool_t PFMET200_HBHE_BeamHaloCleaned_;
  Bool_t PFMETTypeOne200_HBHE_BeamHaloCleaned_;
  Bool_t MET105_IsoTrk50_;
  Bool_t MET120_IsoTrk50_;
  Bool_t SingleJet30_Mu12_SinglePFJet40_;
  Bool_t Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_;
  Bool_t Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_;
  Bool_t Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_;
  Bool_t Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_;
  Bool_t Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t DoublePFJets40_CaloBTagDeepCSV_p71_;
  Bool_t DoublePFJets100_CaloBTagDeepCSV_p71_;
  Bool_t DoublePFJets200_CaloBTagDeepCSV_p71_;
  Bool_t DoublePFJets350_CaloBTagDeepCSV_p71_;
  Bool_t DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t Photon300_NoHE_;
  Bool_t Mu8_TrkIsoVVL_;
  Bool_t Mu8_DiEle12_CaloIdL_TrackIdL_DZ_;
  Bool_t Mu8_DiEle12_CaloIdL_TrackIdL_;
  Bool_t Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_;
  Bool_t Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_;
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_;
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_;
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_;
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_;
  Bool_t Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_;
  Bool_t Mu17_TrkIsoVVL_;
  Bool_t Mu19_TrkIsoVVL_;
  Bool_t BTagMu_AK4DiJet20_Mu5_;
  Bool_t BTagMu_AK4DiJet40_Mu5_;
  Bool_t BTagMu_AK4DiJet70_Mu5_;
  Bool_t BTagMu_AK4DiJet110_Mu5_;
  Bool_t BTagMu_AK4DiJet170_Mu5_;
  Bool_t BTagMu_AK4Jet300_Mu5_;
  Bool_t BTagMu_AK8DiJet170_Mu5_;
  Bool_t BTagMu_AK8Jet170_DoubleMu5_;
  Bool_t BTagMu_AK8Jet300_Mu5_;
  Bool_t BTagMu_AK4DiJet20_Mu5_noalgo_;
  Bool_t BTagMu_AK4DiJet40_Mu5_noalgo_;
  Bool_t BTagMu_AK4DiJet70_Mu5_noalgo_;
  Bool_t BTagMu_AK4DiJet110_Mu5_noalgo_;
  Bool_t BTagMu_AK4DiJet170_Mu5_noalgo_;
  Bool_t BTagMu_AK4Jet300_Mu5_noalgo_;
  Bool_t BTagMu_AK8DiJet170_Mu5_noalgo_;
  Bool_t BTagMu_AK8Jet170_DoubleMu5_noalgo_;
  Bool_t BTagMu_AK8Jet300_Mu5_noalgo_;
  Bool_t Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_;
  Bool_t Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_;
  Bool_t Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t Mu12_DoublePhoton20_;
  Bool_t TriplePhoton_20_20_20_CaloIdLV2_;
  Bool_t TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_;
  Bool_t TriplePhoton_30_30_10_CaloIdLV2_;
  Bool_t TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_;
  Bool_t TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_;
  Bool_t Photon20_;
  Bool_t Photon33_;
  Bool_t Photon50_;
  Bool_t Photon75_;
  Bool_t Photon90_;
  Bool_t Photon120_;
  Bool_t Photon150_;
  Bool_t Photon175_;
  Bool_t Photon200_;
  Bool_t Photon100EB_TightID_TightIso_;
  Bool_t Photon110EB_TightID_TightIso_;
  Bool_t Photon120EB_TightID_TightIso_;
  Bool_t Photon100EBHE10_;
  Bool_t Photon100EEHE10_;
  Bool_t Photon100EE_TightID_TightIso_;
  Bool_t Photon50_R9Id90_HE10_IsoM_;
  Bool_t Photon75_R9Id90_HE10_IsoM_;
  Bool_t Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_;
  Bool_t Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_;
  Bool_t Photon90_R9Id90_HE10_IsoM_;
  Bool_t Photon120_R9Id90_HE10_IsoM_;
  Bool_t Photon165_R9Id90_HE10_IsoM_;
  Bool_t Photon90_CaloIdL_PFHT700_;
  Bool_t Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_;
  Bool_t Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_;
  Bool_t Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_;
  Bool_t Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_;
  Bool_t Photon35_TwoProngs35_;
  Bool_t IsoMu24_TwoProngs35_;
  Bool_t Dimuon0_Jpsi_L1_NoOS_;
  Bool_t Dimuon0_Jpsi_NoVertexing_NoOS_;
  Bool_t Dimuon0_Jpsi_;
  Bool_t Dimuon0_Jpsi_NoVertexing_;
  Bool_t Dimuon0_Jpsi_L1_4R_0er1p5R_;
  Bool_t Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_;
  Bool_t Dimuon0_Jpsi3p5_Muon2_;
  Bool_t Dimuon0_Upsilon_L1_4p5_;
  Bool_t Dimuon0_Upsilon_L1_5_;
  Bool_t Dimuon0_Upsilon_L1_4p5NoOS_;
  Bool_t Dimuon0_Upsilon_L1_4p5er2p0_;
  Bool_t Dimuon0_Upsilon_L1_4p5er2p0M_;
  Bool_t Dimuon0_Upsilon_NoVertexing_;
  Bool_t Dimuon0_Upsilon_L1_5M_;
  Bool_t Dimuon0_LowMass_L1_0er1p5R_;
  Bool_t Dimuon0_LowMass_L1_0er1p5_;
  Bool_t Dimuon0_LowMass_;
  Bool_t Dimuon0_LowMass_L1_4_;
  Bool_t Dimuon0_LowMass_L1_4R_;
  Bool_t Dimuon0_LowMass_L1_TM530_;
  Bool_t Dimuon0_Upsilon_Muon_L1_TM0_;
  Bool_t Dimuon0_Upsilon_Muon_NoL1Mass_;
  Bool_t TripleMu_5_3_3_Mass3p8_DZ_;
  Bool_t TripleMu_10_5_5_DZ_;
  Bool_t TripleMu_12_10_5_;
  Bool_t Tau3Mu_Mu7_Mu1_TkMu1_Tau15_;
  Bool_t Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_;
  Bool_t Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_;
  Bool_t Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_;
  Bool_t DoubleMu3_DZ_PFMET50_PFMHT60_;
  Bool_t DoubleMu3_DZ_PFMET70_PFMHT70_;
  Bool_t DoubleMu3_DZ_PFMET90_PFMHT90_;
  Bool_t DoubleMu3_Trk_Tau3mu_NoL1Mass_;
  Bool_t DoubleMu4_Jpsi_Displaced_;
  Bool_t DoubleMu4_Jpsi_NoVertexing_;
  Bool_t DoubleMu4_JpsiTrkTrk_Displaced_;
  Bool_t DoubleMu43NoFiltersNoVtx_;
  Bool_t DoubleMu48NoFiltersNoVtx_;
  Bool_t Mu43NoFiltersNoVtx_Photon43_CaloIdL_;
  Bool_t Mu48NoFiltersNoVtx_Photon48_CaloIdL_;
  Bool_t Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_;
  Bool_t Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_;
  Bool_t DoubleMu33NoFiltersNoVtxDisplaced_;
  Bool_t DoubleMu40NoFiltersNoVtxDisplaced_;
  Bool_t DoubleMu20_7_Mass0to30_L1_DM4_;
  Bool_t DoubleMu20_7_Mass0to30_L1_DM4EG_;
  Bool_t HT425_;
  Bool_t HT430_DisplacedDijet40_DisplacedTrack_;
  Bool_t HT500_DisplacedDijet40_DisplacedTrack_;
  Bool_t HT430_DisplacedDijet60_DisplacedTrack_;
  Bool_t HT400_DisplacedDijet40_DisplacedTrack_;
  Bool_t HT650_DisplacedDijet60_Inclusive_;
  Bool_t HT550_DisplacedDijet60_Inclusive_;
  Bool_t DiJet110_35_Mjj650_PFMET110_;
  Bool_t DiJet110_35_Mjj650_PFMET120_;
  Bool_t DiJet110_35_Mjj650_PFMET130_;
  Bool_t TripleJet110_35_35_Mjj650_PFMET110_;
  Bool_t TripleJet110_35_35_Mjj650_PFMET120_;
  Bool_t TripleJet110_35_35_Mjj650_PFMET130_;
  Bool_t Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_;
  Bool_t Ele28_eta2p1_WPTight_Gsf_HT150_;
  Bool_t Ele28_HighEta_SC20_Mass55_;
  Bool_t DoubleMu20_7_Mass0to30_Photon23_;
  Bool_t Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_;
  Bool_t Ele15_IsoVVVL_PFHT450_PFMET50_;
  Bool_t Ele15_IsoVVVL_PFHT450_;
  Bool_t Ele50_IsoVVVL_PFHT450_;
  Bool_t Ele15_IsoVVVL_PFHT600_;
  Bool_t Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_;
  Bool_t Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_;
  Bool_t Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_;
  Bool_t Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_;
  Bool_t Mu15_IsoVVVL_PFHT450_PFMET50_;
  Bool_t Mu15_IsoVVVL_PFHT450_;
  Bool_t Mu50_IsoVVVL_PFHT450_;
  Bool_t Mu15_IsoVVVL_PFHT600_;
  Bool_t Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_;
  Bool_t Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_;
  Bool_t Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_;
  Bool_t Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_;
  Bool_t Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_;
  Bool_t Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_;
  Bool_t Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_;
  Bool_t Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_;
  Bool_t Dimuon10_PsiPrime_Barrel_Seagulls_;
  Bool_t Dimuon20_Jpsi_Barrel_Seagulls_;
  Bool_t Dimuon12_Upsilon_y1p4_;
  Bool_t Dimuon14_Phi_Barrel_Seagulls_;
  Bool_t Dimuon18_PsiPrime_;
  Bool_t Dimuon25_Jpsi_;
  Bool_t Dimuon18_PsiPrime_noCorrL1_;
  Bool_t Dimuon24_Upsilon_noCorrL1_;
  Bool_t Dimuon24_Phi_noCorrL1_;
  Bool_t Dimuon25_Jpsi_noCorrL1_;
  Bool_t DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_;
  Bool_t DiMu9_Ele9_CaloIdL_TrackIdL_DZ_;
  Bool_t DiMu9_Ele9_CaloIdL_TrackIdL_;
  Bool_t DoubleIsoMu20_eta2p1_;
  Bool_t TrkMu12_DoubleTrkMu5NoFiltersNoVtx_;
  Bool_t TrkMu16_DoubleTrkMu6NoFiltersNoVtx_;
  Bool_t TrkMu17_DoubleTrkMu8NoFiltersNoVtx_;
  Bool_t Mu8_;
  Bool_t Mu17_;
  Bool_t Mu19_;
  Bool_t Mu17_Photon30_IsoCaloId_;
  Bool_t Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t Ele8_CaloIdM_TrackIdM_PFJet30_;
  Bool_t Ele17_CaloIdM_TrackIdM_PFJet30_;
  Bool_t Ele23_CaloIdM_TrackIdM_PFJet30_;
  Bool_t Ele50_CaloIdVT_GsfTrkIdT_PFJet165_;
  Bool_t Ele115_CaloIdVT_GsfTrkIdT_;
  Bool_t Ele135_CaloIdVT_GsfTrkIdT_;
  Bool_t Ele145_CaloIdVT_GsfTrkIdT_;
  Bool_t Ele200_CaloIdVT_GsfTrkIdT_;
  Bool_t Ele250_CaloIdVT_GsfTrkIdT_;
  Bool_t Ele300_CaloIdVT_GsfTrkIdT_;
  Bool_t PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_;
  Bool_t PFHT330PT30_QuadPFJet_75_60_45_40_;
  Bool_t PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_;
  Bool_t PFHT400_SixPFJet32_;
  Bool_t PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_;
  Bool_t PFHT450_SixPFJet36_;
  Bool_t PFHT350_;
  Bool_t PFHT350MinPFJet15_;
  Bool_t Photon60_R9Id90_CaloIdL_IsoL_;
  Bool_t Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_;
  Bool_t Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_;
  Bool_t ECALHT800_;
  Bool_t DiSC30_18_EIso_AND_HE_Mass70_;
  Bool_t Physics_;
  Bool_t Physics_part0_;
  Bool_t Physics_part1_;
  Bool_t Physics_part2_;
  Bool_t Physics_part3_;
  Bool_t Physics_part4_;
  Bool_t Physics_part5_;
  Bool_t Physics_part6_;
  Bool_t Physics_part7_;
  Bool_t Random_;
  Bool_t ZeroBias_;
  Bool_t ZeroBias_Alignment_;
  Bool_t ZeroBias_part0_;
  Bool_t ZeroBias_part1_;
  Bool_t ZeroBias_part2_;
  Bool_t ZeroBias_part3_;
  Bool_t ZeroBias_part4_;
  Bool_t ZeroBias_part5_;
  Bool_t ZeroBias_part6_;
  Bool_t ZeroBias_part7_;
  Bool_t AK4CaloJet30_;
  Bool_t AK4CaloJet40_;
  Bool_t AK4CaloJet50_;
  Bool_t AK4CaloJet80_;
  Bool_t AK4CaloJet100_;
  Bool_t AK4CaloJet120_;
  Bool_t AK4PFJet30_;
  Bool_t AK4PFJet50_;
  Bool_t AK4PFJet80_;
  Bool_t AK4PFJet100_;
  Bool_t AK4PFJet120_;
  Bool_t SinglePhoton10_Eta3p1ForPPRef_;
  Bool_t SinglePhoton20_Eta3p1ForPPRef_;
  Bool_t SinglePhoton30_Eta3p1ForPPRef_;
  Bool_t Photon20_HoverELoose_;
  Bool_t Photon30_HoverELoose_;
  Bool_t EcalCalibration_;
  Bool_t HcalCalibration_;
  Bool_t L1UnpairedBunchBptxMinus_;
  Bool_t L1UnpairedBunchBptxPlus_;
  Bool_t L1NotBptxOR_;
  Bool_t L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_;
  Bool_t CDC_L2cosmic_5_er1p0_;
  Bool_t CDC_L2cosmic_5p5_er1p0_;
  Bool_t HcalNZS_;
  Bool_t HcalPhiSym_;
  Bool_t HcalIsolatedbunch_;
  Bool_t IsoTrackHB_;
  Bool_t IsoTrackHE_;
  Bool_t ZeroBias_FirstCollisionAfterAbortGap_;
  Bool_t ZeroBias_IsolatedBunches_;
  Bool_t ZeroBias_FirstCollisionInTrain_;
  Bool_t ZeroBias_LastCollisionInTrain_;
  Bool_t ZeroBias_FirstBXAfterTrain_;
  Bool_t IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_;
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_;
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_;
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_;
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_;
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_;
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_;
  Bool_t MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_;
  Bool_t MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_;
  Bool_t MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_;
  Bool_t MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_;
  Bool_t MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_;
  Bool_t Ele16_Ele12_Ele8_CaloIdL_TrackIdL_;
  Bool_t Rsq0p35_;
  Bool_t Rsq0p40_;
  Bool_t RsqMR300_Rsq0p09_MR200_;
  Bool_t RsqMR320_Rsq0p09_MR200_;
  Bool_t RsqMR300_Rsq0p09_MR200_4jet_;
  Bool_t RsqMR320_Rsq0p09_MR200_4jet_;
  Bool_t IsoMu27_MET90_;
  Bool_t DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_;
  Bool_t DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_;
  Bool_t DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_;
  Bool_t DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_;
  Bool_t DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_;
  Bool_t DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_;
  Bool_t DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_;
  Bool_t DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_;
  Bool_t VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_;
  Bool_t VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_;
  Bool_t VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_;
  Bool_t Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_;
  Bool_t Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_;
  Bool_t Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_;
  Bool_t PFMET100_PFMHT100_IDTight_PFHT60_;
  Bool_t PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_;
  Bool_t PFMETTypeOne100_PFMHT100_IDTight_PFHT60_;
  Bool_t Mu18_Mu9_SameSign_;
  Bool_t Mu18_Mu9_SameSign_DZ_;
  Bool_t Mu18_Mu9_;
  Bool_t Mu18_Mu9_DZ_;
  Bool_t Mu20_Mu10_SameSign_;
  Bool_t Mu20_Mu10_SameSign_DZ_;
  Bool_t Mu20_Mu10_;
  Bool_t Mu20_Mu10_DZ_;
  Bool_t Mu23_Mu12_SameSign_;
  Bool_t Mu23_Mu12_SameSign_DZ_;
  Bool_t Mu23_Mu12_;
  Bool_t Mu23_Mu12_DZ_;
  Bool_t DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_;
  Bool_t DoubleMu2_Jpsi_DoubleTkMu0_Phi_;
  Bool_t DoubleMu3_DCA_PFMET50_PFMHT60_;
  Bool_t TripleMu_5_3_3_Mass3p8_DCA_;
  Bool_t QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;
  Bool_t QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;
  Bool_t QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;
  Bool_t QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_;
  Bool_t QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_;
  Bool_t QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_;
  Bool_t QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_;
  Bool_t QuadPFJet98_83_71_15_;
  Bool_t QuadPFJet103_88_75_15_;
  Bool_t QuadPFJet105_88_76_15_;
  Bool_t QuadPFJet111_90_80_15_;
  Bool_t AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_;
  Bool_t AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_;
  Bool_t AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_;
  Bool_t AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_;
  Bool_t AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_;
  Bool_t Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_;
  Bool_t Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_;
  Bool_t Mu12_IP6_part0_;
  Bool_t Mu12_IP6_part1_;
  Bool_t Mu12_IP6_part2_;
  Bool_t Mu12_IP6_part3_;
  Bool_t Mu12_IP6_part4_;
  Bool_t Mu9_IP5_part0_;
  Bool_t Mu9_IP5_part1_;
  Bool_t Mu9_IP5_part2_;
  Bool_t Mu9_IP5_part3_;
  Bool_t Mu9_IP5_part4_;
  Bool_t Mu7_IP4_part0_;
  Bool_t Mu7_IP4_part1_;
  Bool_t Mu7_IP4_part2_;
  Bool_t Mu7_IP4_part3_;
  Bool_t Mu7_IP4_part4_;
  Bool_t Mu9_IP4_part0_;
  Bool_t Mu9_IP4_part1_;
  Bool_t Mu9_IP4_part2_;
  Bool_t Mu9_IP4_part3_;
  Bool_t Mu9_IP4_part4_;
  Bool_t Mu8_IP5_part0_;
  Bool_t Mu8_IP5_part1_;
  Bool_t Mu8_IP5_part2_;
  Bool_t Mu8_IP5_part3_;
  Bool_t Mu8_IP5_part4_;
  Bool_t Mu8_IP6_part0_;
  Bool_t Mu8_IP6_part1_;
  Bool_t Mu8_IP6_part2_;
  Bool_t Mu8_IP6_part3_;
  Bool_t Mu8_IP6_part4_;
  Bool_t Mu9_IP6_part0_;
  Bool_t Mu9_IP6_part1_;
  Bool_t Mu9_IP6_part2_;
  Bool_t Mu9_IP6_part3_;
  Bool_t Mu9_IP6_part4_;
  Bool_t Mu8_IP3_part0_;
  Bool_t Mu8_IP3_part1_;
  Bool_t Mu8_IP3_part2_;
  Bool_t Mu8_IP3_part3_;
  Bool_t Mu8_IP3_part4_;
  Bool_t QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;
  Bool_t TrkMu6NoFiltersNoVtx_;
  Bool_t TrkMu16NoFiltersNoVtx_;
  Bool_t DoubleTrkMu_16_6_NoFiltersNoVtx_;
  Bool_t HLTriggerFinalPath_;
  void setHLTriggerFirstPath(const Bool_t value) {HLTriggerFirstPath_ = value;}
  void setAK8PFJet360_TrimMass30(const Bool_t value) {AK8PFJet360_TrimMass30_ = value;}
  void setAK8PFJet380_TrimMass30(const Bool_t value) {AK8PFJet380_TrimMass30_ = value;}
  void setAK8PFJet400_TrimMass30(const Bool_t value) {AK8PFJet400_TrimMass30_ = value;}
  void setAK8PFJet420_TrimMass30(const Bool_t value) {AK8PFJet420_TrimMass30_ = value;}
  void setAK8PFHT750_TrimMass50(const Bool_t value) {AK8PFHT750_TrimMass50_ = value;}
  void setAK8PFHT800_TrimMass50(const Bool_t value) {AK8PFHT800_TrimMass50_ = value;}
  void setAK8PFHT850_TrimMass50(const Bool_t value) {AK8PFHT850_TrimMass50_ = value;}
  void setAK8PFHT900_TrimMass50(const Bool_t value) {AK8PFHT900_TrimMass50_ = value;}
  void setCaloJet500_NoJetID(const Bool_t value) {CaloJet500_NoJetID_ = value;}
  void setCaloJet550_NoJetID(const Bool_t value) {CaloJet550_NoJetID_ = value;}
  void setDoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL(const Bool_t value) {DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_ = value;}
  void setDoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon(const Bool_t value) {DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_ = value;}
  void setTrimuon5_3p5_2_Upsilon_Muon(const Bool_t value) {Trimuon5_3p5_2_Upsilon_Muon_ = value;}
  void setTrimuonOpen_5_3p5_2_Upsilon_Muon(const Bool_t value) {TrimuonOpen_5_3p5_2_Upsilon_Muon_ = value;}
  void setDoubleEle25_CaloIdL_MW(const Bool_t value) {DoubleEle25_CaloIdL_MW_ = value;}
  void setDoubleEle27_CaloIdL_MW(const Bool_t value) {DoubleEle27_CaloIdL_MW_ = value;}
  void setDoubleEle33_CaloIdL_MW(const Bool_t value) {DoubleEle33_CaloIdL_MW_ = value;}
  void setDoubleEle24_eta2p1_WPTight_Gsf(const Bool_t value) {DoubleEle24_eta2p1_WPTight_Gsf_ = value;}
  void setDoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350(const Bool_t value) {DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_ = value;}
  void setDoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350(const Bool_t value) {DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_ = value;}
  void setEle27_Ele37_CaloIdL_MW(const Bool_t value) {Ele27_Ele37_CaloIdL_MW_ = value;}
  void setMu27_Ele37_CaloIdL_MW(const Bool_t value) {Mu27_Ele37_CaloIdL_MW_ = value;}
  void setMu37_Ele27_CaloIdL_MW(const Bool_t value) {Mu37_Ele27_CaloIdL_MW_ = value;}
  void setMu37_TkMu27(const Bool_t value) {Mu37_TkMu27_ = value;}
  void setDoubleMu4_3_Bs(const Bool_t value) {DoubleMu4_3_Bs_ = value;}
  void setDoubleMu4_3_Jpsi(const Bool_t value) {DoubleMu4_3_Jpsi_ = value;}
  void setDoubleMu4_JpsiTrk_Displaced(const Bool_t value) {DoubleMu4_JpsiTrk_Displaced_ = value;}
  void setDoubleMu4_LowMassNonResonantTrk_Displaced(const Bool_t value) {DoubleMu4_LowMassNonResonantTrk_Displaced_ = value;}
  void setDoubleMu3_Trk_Tau3mu(const Bool_t value) {DoubleMu3_Trk_Tau3mu_ = value;}
  void setDoubleMu3_TkMu_DsTau3Mu(const Bool_t value) {DoubleMu3_TkMu_DsTau3Mu_ = value;}
  void setDoubleMu4_PsiPrimeTrk_Displaced(const Bool_t value) {DoubleMu4_PsiPrimeTrk_Displaced_ = value;}
  void setDoubleMu4_Mass3p8_DZ_PFHT350(const Bool_t value) {DoubleMu4_Mass3p8_DZ_PFHT350_ = value;}
  void setMu3_PFJet40(const Bool_t value) {Mu3_PFJet40_ = value;}
  void setMu7p5_L2Mu2_Jpsi(const Bool_t value) {Mu7p5_L2Mu2_Jpsi_ = value;}
  void setMu7p5_L2Mu2_Upsilon(const Bool_t value) {Mu7p5_L2Mu2_Upsilon_ = value;}
  void setMu7p5_Track2_Jpsi(const Bool_t value) {Mu7p5_Track2_Jpsi_ = value;}
  void setMu7p5_Track3p5_Jpsi(const Bool_t value) {Mu7p5_Track3p5_Jpsi_ = value;}
  void setMu7p5_Track7_Jpsi(const Bool_t value) {Mu7p5_Track7_Jpsi_ = value;}
  void setMu7p5_Track2_Upsilon(const Bool_t value) {Mu7p5_Track2_Upsilon_ = value;}
  void setMu7p5_Track3p5_Upsilon(const Bool_t value) {Mu7p5_Track3p5_Upsilon_ = value;}
  void setMu7p5_Track7_Upsilon(const Bool_t value) {Mu7p5_Track7_Upsilon_ = value;}
  void setMu3_L1SingleMu5orSingleMu7(const Bool_t value) {Mu3_L1SingleMu5orSingleMu7_ = value;}
  void setDoublePhoton33_CaloIdL(const Bool_t value) {DoublePhoton33_CaloIdL_ = value;}
  void setDoublePhoton70(const Bool_t value) {DoublePhoton70_ = value;}
  void setDoublePhoton85(const Bool_t value) {DoublePhoton85_ = value;}
  void setEle20_WPTight_Gsf(const Bool_t value) {Ele20_WPTight_Gsf_ = value;}
  void setEle15_WPLoose_Gsf(const Bool_t value) {Ele15_WPLoose_Gsf_ = value;}
  void setEle17_WPLoose_Gsf(const Bool_t value) {Ele17_WPLoose_Gsf_ = value;}
  void setEle20_WPLoose_Gsf(const Bool_t value) {Ele20_WPLoose_Gsf_ = value;}
  void setEle20_eta2p1_WPLoose_Gsf(const Bool_t value) {Ele20_eta2p1_WPLoose_Gsf_ = value;}
  void setDiEle27_WPTightCaloOnly_L1DoubleEG(const Bool_t value) {DiEle27_WPTightCaloOnly_L1DoubleEG_ = value;}
  void setEle27_WPTight_Gsf(const Bool_t value) {Ele27_WPTight_Gsf_ = value;}
  void setEle28_WPTight_Gsf(const Bool_t value) {Ele28_WPTight_Gsf_ = value;}
  void setEle30_WPTight_Gsf(const Bool_t value) {Ele30_WPTight_Gsf_ = value;}
  void setEle32_WPTight_Gsf(const Bool_t value) {Ele32_WPTight_Gsf_ = value;}
  void setEle35_WPTight_Gsf(const Bool_t value) {Ele35_WPTight_Gsf_ = value;}
  void setEle35_WPTight_Gsf_L1EGMT(const Bool_t value) {Ele35_WPTight_Gsf_L1EGMT_ = value;}
  void setEle38_WPTight_Gsf(const Bool_t value) {Ele38_WPTight_Gsf_ = value;}
  void setEle40_WPTight_Gsf(const Bool_t value) {Ele40_WPTight_Gsf_ = value;}
  void setEle32_WPTight_Gsf_L1DoubleEG(const Bool_t value) {Ele32_WPTight_Gsf_L1DoubleEG_ = value;}
  void setEle24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1(const Bool_t value) {Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_ = value;}
  void setEle24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1(const Bool_t value) {Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_ = value;}
  void setEle24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1(const Bool_t value) {Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_ = value;}
  void setEle24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1(const Bool_t value) {Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_ = value;}
  void setEle24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1(const Bool_t value) {Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_ = value;}
  void setEle24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1(const Bool_t value) {Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_ = value;}
  void setHT450_Beamspot(const Bool_t value) {HT450_Beamspot_ = value;}
  void setHT300_Beamspot(const Bool_t value) {HT300_Beamspot_ = value;}
  void setZeroBias_Beamspot(const Bool_t value) {ZeroBias_Beamspot_ = value;}
  void setIsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1(const Bool_t value) {IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_ = value;}
  void setIsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1(const Bool_t value) {IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_ = value;}
  void setIsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1(const Bool_t value) {IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_ = value;}
  void setIsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1(const Bool_t value) {IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_ = value;}
  void setIsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1(const Bool_t value) {IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_ = value;}
  void setIsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1(const Bool_t value) {IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_ = value;}
  void setIsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1(const Bool_t value) {IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_ = value;}
  void setIsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1(const Bool_t value) {IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_ = value;}
  void setIsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1(const Bool_t value) {IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_ = value;}
  void setIsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1(const Bool_t value) {IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_ = value;}
  void setIsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1(const Bool_t value) {IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_ = value;}
  void setIsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1(const Bool_t value) {IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_ = value;}
  void setIsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1(const Bool_t value) {IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_ = value;}
  void setIsoMu20(const Bool_t value) {IsoMu20_ = value;}
  void setIsoMu24(const Bool_t value) {IsoMu24_ = value;}
  void setIsoMu24_eta2p1(const Bool_t value) {IsoMu24_eta2p1_ = value;}
  void setIsoMu27(const Bool_t value) {IsoMu27_ = value;}
  void setIsoMu30(const Bool_t value) {IsoMu30_ = value;}
  void setUncorrectedJetE30_NoBPTX(const Bool_t value) {UncorrectedJetE30_NoBPTX_ = value;}
  void setUncorrectedJetE30_NoBPTX3BX(const Bool_t value) {UncorrectedJetE30_NoBPTX3BX_ = value;}
  void setUncorrectedJetE60_NoBPTX3BX(const Bool_t value) {UncorrectedJetE60_NoBPTX3BX_ = value;}
  void setUncorrectedJetE70_NoBPTX3BX(const Bool_t value) {UncorrectedJetE70_NoBPTX3BX_ = value;}
  void setL1SingleMu18(const Bool_t value) {L1SingleMu18_ = value;}
  void setL1SingleMu25(const Bool_t value) {L1SingleMu25_ = value;}
  void setL2Mu10(const Bool_t value) {L2Mu10_ = value;}
  void setL2Mu10_NoVertex_NoBPTX3BX(const Bool_t value) {L2Mu10_NoVertex_NoBPTX3BX_ = value;}
  void setL2Mu10_NoVertex_NoBPTX(const Bool_t value) {L2Mu10_NoVertex_NoBPTX_ = value;}
  void setL2Mu45_NoVertex_3Sta_NoBPTX3BX(const Bool_t value) {L2Mu45_NoVertex_3Sta_NoBPTX3BX_ = value;}
  void setL2Mu40_NoVertex_3Sta_NoBPTX3BX(const Bool_t value) {L2Mu40_NoVertex_3Sta_NoBPTX3BX_ = value;}
  void setL2Mu50(const Bool_t value) {L2Mu50_ = value;}
  void setL2Mu23NoVtx_2Cha(const Bool_t value) {L2Mu23NoVtx_2Cha_ = value;}
  void setL2Mu23NoVtx_2Cha_CosmicSeed(const Bool_t value) {L2Mu23NoVtx_2Cha_CosmicSeed_ = value;}
  void setDoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4(const Bool_t value) {DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_ = value;}
  void setDoubleL2Mu30NoVtx_2Cha_Eta2p4(const Bool_t value) {DoubleL2Mu30NoVtx_2Cha_Eta2p4_ = value;}
  void setDoubleL2Mu50(const Bool_t value) {DoubleL2Mu50_ = value;}
  void setDoubleL2Mu23NoVtx_2Cha_CosmicSeed(const Bool_t value) {DoubleL2Mu23NoVtx_2Cha_CosmicSeed_ = value;}
  void setDoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched(const Bool_t value) {DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_ = value;}
  void setDoubleL2Mu25NoVtx_2Cha_CosmicSeed(const Bool_t value) {DoubleL2Mu25NoVtx_2Cha_CosmicSeed_ = value;}
  void setDoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched(const Bool_t value) {DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_ = value;}
  void setDoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4(const Bool_t value) {DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_ = value;}
  void setDoubleL2Mu23NoVtx_2Cha(const Bool_t value) {DoubleL2Mu23NoVtx_2Cha_ = value;}
  void setDoubleL2Mu23NoVtx_2Cha_NoL2Matched(const Bool_t value) {DoubleL2Mu23NoVtx_2Cha_NoL2Matched_ = value;}
  void setDoubleL2Mu25NoVtx_2Cha(const Bool_t value) {DoubleL2Mu25NoVtx_2Cha_ = value;}
  void setDoubleL2Mu25NoVtx_2Cha_NoL2Matched(const Bool_t value) {DoubleL2Mu25NoVtx_2Cha_NoL2Matched_ = value;}
  void setDoubleL2Mu25NoVtx_2Cha_Eta2p4(const Bool_t value) {DoubleL2Mu25NoVtx_2Cha_Eta2p4_ = value;}
  void setMu17_TrkIsoVVL_Mu8_TrkIsoVVL(const Bool_t value) {Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_ = value;}
  void setMu19_TrkIsoVVL_Mu9_TrkIsoVVL(const Bool_t value) {Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_ = value;}
  void setMu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ(const Bool_t value) {Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_ = value;}
  void setMu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ(const Bool_t value) {Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_ = value;}
  void setMu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8(const Bool_t value) {Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_ = value;}
  void setMu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8(const Bool_t value) {Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_ = value;}
  void setMu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8(const Bool_t value) {Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_ = value;}
  void setMu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8(const Bool_t value) {Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_ = value;}
  void setMu25_TkMu0_Onia(const Bool_t value) {Mu25_TkMu0_Onia_ = value;}
  void setMu30_TkMu0_Psi(const Bool_t value) {Mu30_TkMu0_Psi_ = value;}
  void setMu30_TkMu0_Upsilon(const Bool_t value) {Mu30_TkMu0_Upsilon_ = value;}
  void setMu20_TkMu0_Phi(const Bool_t value) {Mu20_TkMu0_Phi_ = value;}
  void setMu25_TkMu0_Phi(const Bool_t value) {Mu25_TkMu0_Phi_ = value;}
  void setMu12(const Bool_t value) {Mu12_ = value;}
  void setMu15(const Bool_t value) {Mu15_ = value;}
  void setMu20(const Bool_t value) {Mu20_ = value;}
  void setMu27(const Bool_t value) {Mu27_ = value;}
  void setMu50(const Bool_t value) {Mu50_ = value;}
  void setMu55(const Bool_t value) {Mu55_ = value;}
  void setOldMu100(const Bool_t value) {OldMu100_ = value;}
  void setTkMu100(const Bool_t value) {TkMu100_ = value;}
  void setDiPFJetAve40(const Bool_t value) {DiPFJetAve40_ = value;}
  void setDiPFJetAve60(const Bool_t value) {DiPFJetAve60_ = value;}
  void setDiPFJetAve80(const Bool_t value) {DiPFJetAve80_ = value;}
  void setDiPFJetAve140(const Bool_t value) {DiPFJetAve140_ = value;}
  void setDiPFJetAve200(const Bool_t value) {DiPFJetAve200_ = value;}
  void setDiPFJetAve260(const Bool_t value) {DiPFJetAve260_ = value;}
  void setDiPFJetAve320(const Bool_t value) {DiPFJetAve320_ = value;}
  void setDiPFJetAve400(const Bool_t value) {DiPFJetAve400_ = value;}
  void setDiPFJetAve500(const Bool_t value) {DiPFJetAve500_ = value;}
  void setDiPFJetAve60_HFJEC(const Bool_t value) {DiPFJetAve60_HFJEC_ = value;}
  void setDiPFJetAve80_HFJEC(const Bool_t value) {DiPFJetAve80_HFJEC_ = value;}
  void setDiPFJetAve100_HFJEC(const Bool_t value) {DiPFJetAve100_HFJEC_ = value;}
  void setDiPFJetAve160_HFJEC(const Bool_t value) {DiPFJetAve160_HFJEC_ = value;}
  void setDiPFJetAve220_HFJEC(const Bool_t value) {DiPFJetAve220_HFJEC_ = value;}
  void setDiPFJetAve300_HFJEC(const Bool_t value) {DiPFJetAve300_HFJEC_ = value;}
  void setAK8PFJet15(const Bool_t value) {AK8PFJet15_ = value;}
  void setAK8PFJet25(const Bool_t value) {AK8PFJet25_ = value;}
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
  void setAK8PFJet550(const Bool_t value) {AK8PFJet550_ = value;}
  void setPFJet15(const Bool_t value) {PFJet15_ = value;}
  void setPFJet25(const Bool_t value) {PFJet25_ = value;}
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
  void setPFJet550(const Bool_t value) {PFJet550_ = value;}
  void setPFJetFwd15(const Bool_t value) {PFJetFwd15_ = value;}
  void setPFJetFwd25(const Bool_t value) {PFJetFwd25_ = value;}
  void setPFJetFwd40(const Bool_t value) {PFJetFwd40_ = value;}
  void setPFJetFwd60(const Bool_t value) {PFJetFwd60_ = value;}
  void setPFJetFwd80(const Bool_t value) {PFJetFwd80_ = value;}
  void setPFJetFwd140(const Bool_t value) {PFJetFwd140_ = value;}
  void setPFJetFwd200(const Bool_t value) {PFJetFwd200_ = value;}
  void setPFJetFwd260(const Bool_t value) {PFJetFwd260_ = value;}
  void setPFJetFwd320(const Bool_t value) {PFJetFwd320_ = value;}
  void setPFJetFwd400(const Bool_t value) {PFJetFwd400_ = value;}
  void setPFJetFwd450(const Bool_t value) {PFJetFwd450_ = value;}
  void setPFJetFwd500(const Bool_t value) {PFJetFwd500_ = value;}
  void setAK8PFJetFwd15(const Bool_t value) {AK8PFJetFwd15_ = value;}
  void setAK8PFJetFwd25(const Bool_t value) {AK8PFJetFwd25_ = value;}
  void setAK8PFJetFwd40(const Bool_t value) {AK8PFJetFwd40_ = value;}
  void setAK8PFJetFwd60(const Bool_t value) {AK8PFJetFwd60_ = value;}
  void setAK8PFJetFwd80(const Bool_t value) {AK8PFJetFwd80_ = value;}
  void setAK8PFJetFwd140(const Bool_t value) {AK8PFJetFwd140_ = value;}
  void setAK8PFJetFwd200(const Bool_t value) {AK8PFJetFwd200_ = value;}
  void setAK8PFJetFwd260(const Bool_t value) {AK8PFJetFwd260_ = value;}
  void setAK8PFJetFwd320(const Bool_t value) {AK8PFJetFwd320_ = value;}
  void setAK8PFJetFwd400(const Bool_t value) {AK8PFJetFwd400_ = value;}
  void setAK8PFJetFwd450(const Bool_t value) {AK8PFJetFwd450_ = value;}
  void setAK8PFJetFwd500(const Bool_t value) {AK8PFJetFwd500_ = value;}
  void setPFHT180(const Bool_t value) {PFHT180_ = value;}
  void setPFHT250(const Bool_t value) {PFHT250_ = value;}
  void setPFHT370(const Bool_t value) {PFHT370_ = value;}
  void setPFHT430(const Bool_t value) {PFHT430_ = value;}
  void setPFHT510(const Bool_t value) {PFHT510_ = value;}
  void setPFHT590(const Bool_t value) {PFHT590_ = value;}
  void setPFHT680(const Bool_t value) {PFHT680_ = value;}
  void setPFHT780(const Bool_t value) {PFHT780_ = value;}
  void setPFHT890(const Bool_t value) {PFHT890_ = value;}
  void setPFHT1050(const Bool_t value) {PFHT1050_ = value;}
  void setPFHT500_PFMET100_PFMHT100_IDTight(const Bool_t value) {PFHT500_PFMET100_PFMHT100_IDTight_ = value;}
  void setPFHT500_PFMET110_PFMHT110_IDTight(const Bool_t value) {PFHT500_PFMET110_PFMHT110_IDTight_ = value;}
  void setPFHT700_PFMET85_PFMHT85_IDTight(const Bool_t value) {PFHT700_PFMET85_PFMHT85_IDTight_ = value;}
  void setPFHT700_PFMET95_PFMHT95_IDTight(const Bool_t value) {PFHT700_PFMET95_PFMHT95_IDTight_ = value;}
  void setPFHT800_PFMET75_PFMHT75_IDTight(const Bool_t value) {PFHT800_PFMET75_PFMHT75_IDTight_ = value;}
  void setPFHT800_PFMET85_PFMHT85_IDTight(const Bool_t value) {PFHT800_PFMET85_PFMHT85_IDTight_ = value;}
  void setPFMET110_PFMHT110_IDTight(const Bool_t value) {PFMET110_PFMHT110_IDTight_ = value;}
  void setPFMET120_PFMHT120_IDTight(const Bool_t value) {PFMET120_PFMHT120_IDTight_ = value;}
  void setPFMET130_PFMHT130_IDTight(const Bool_t value) {PFMET130_PFMHT130_IDTight_ = value;}
  void setPFMET140_PFMHT140_IDTight(const Bool_t value) {PFMET140_PFMHT140_IDTight_ = value;}
  void setPFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1(const Bool_t value) {PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_ = value;}
  void setPFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1(const Bool_t value) {PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_ = value;}
  void setPFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1(const Bool_t value) {PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_ = value;}
  void setPFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1(const Bool_t value) {PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_ = value;}
  void setPFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1(const Bool_t value) {PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_ = value;}
  void setPFMET120_PFMHT120_IDTight_PFHT60(const Bool_t value) {PFMET120_PFMHT120_IDTight_PFHT60_ = value;}
  void setPFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60(const Bool_t value) {PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ = value;}
  void setPFMETTypeOne120_PFMHT120_IDTight_PFHT60(const Bool_t value) {PFMETTypeOne120_PFMHT120_IDTight_PFHT60_ = value;}
  void setPFMETTypeOne110_PFMHT110_IDTight(const Bool_t value) {PFMETTypeOne110_PFMHT110_IDTight_ = value;}
  void setPFMETTypeOne120_PFMHT120_IDTight(const Bool_t value) {PFMETTypeOne120_PFMHT120_IDTight_ = value;}
  void setPFMETTypeOne130_PFMHT130_IDTight(const Bool_t value) {PFMETTypeOne130_PFMHT130_IDTight_ = value;}
  void setPFMETTypeOne140_PFMHT140_IDTight(const Bool_t value) {PFMETTypeOne140_PFMHT140_IDTight_ = value;}
  void setPFMETNoMu110_PFMHTNoMu110_IDTight(const Bool_t value) {PFMETNoMu110_PFMHTNoMu110_IDTight_ = value;}
  void setPFMETNoMu120_PFMHTNoMu120_IDTight(const Bool_t value) {PFMETNoMu120_PFMHTNoMu120_IDTight_ = value;}
  void setPFMETNoMu130_PFMHTNoMu130_IDTight(const Bool_t value) {PFMETNoMu130_PFMHTNoMu130_IDTight_ = value;}
  void setPFMETNoMu140_PFMHTNoMu140_IDTight(const Bool_t value) {PFMETNoMu140_PFMHTNoMu140_IDTight_ = value;}
  void setMonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight(const Bool_t value) {MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_ = value;}
  void setMonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight(const Bool_t value) {MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_ = value;}
  void setMonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight(const Bool_t value) {MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_ = value;}
  void setMonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight(const Bool_t value) {MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_ = value;}
  void setL1ETMHadSeeds(const Bool_t value) {L1ETMHadSeeds_ = value;}
  void setCaloMHT90(const Bool_t value) {CaloMHT90_ = value;}
  void setCaloMET80_NotCleaned(const Bool_t value) {CaloMET80_NotCleaned_ = value;}
  void setCaloMET90_NotCleaned(const Bool_t value) {CaloMET90_NotCleaned_ = value;}
  void setCaloMET100_NotCleaned(const Bool_t value) {CaloMET100_NotCleaned_ = value;}
  void setCaloMET110_NotCleaned(const Bool_t value) {CaloMET110_NotCleaned_ = value;}
  void setCaloMET250_NotCleaned(const Bool_t value) {CaloMET250_NotCleaned_ = value;}
  void setCaloMET70_HBHECleaned(const Bool_t value) {CaloMET70_HBHECleaned_ = value;}
  void setCaloMET80_HBHECleaned(const Bool_t value) {CaloMET80_HBHECleaned_ = value;}
  void setCaloMET90_HBHECleaned(const Bool_t value) {CaloMET90_HBHECleaned_ = value;}
  void setCaloMET100_HBHECleaned(const Bool_t value) {CaloMET100_HBHECleaned_ = value;}
  void setCaloMET250_HBHECleaned(const Bool_t value) {CaloMET250_HBHECleaned_ = value;}
  void setCaloMET300_HBHECleaned(const Bool_t value) {CaloMET300_HBHECleaned_ = value;}
  void setCaloMET350_HBHECleaned(const Bool_t value) {CaloMET350_HBHECleaned_ = value;}
  void setPFMET200_NotCleaned(const Bool_t value) {PFMET200_NotCleaned_ = value;}
  void setPFMET200_HBHECleaned(const Bool_t value) {PFMET200_HBHECleaned_ = value;}
  void setPFMET250_HBHECleaned(const Bool_t value) {PFMET250_HBHECleaned_ = value;}
  void setPFMET300_HBHECleaned(const Bool_t value) {PFMET300_HBHECleaned_ = value;}
  void setPFMET200_HBHE_BeamHaloCleaned(const Bool_t value) {PFMET200_HBHE_BeamHaloCleaned_ = value;}
  void setPFMETTypeOne200_HBHE_BeamHaloCleaned(const Bool_t value) {PFMETTypeOne200_HBHE_BeamHaloCleaned_ = value;}
  void setMET105_IsoTrk50(const Bool_t value) {MET105_IsoTrk50_ = value;}
  void setMET120_IsoTrk50(const Bool_t value) {MET120_IsoTrk50_ = value;}
  void setSingleJet30_Mu12_SinglePFJet40(const Bool_t value) {SingleJet30_Mu12_SinglePFJet40_ = value;}
  void setMu12_DoublePFJets40_CaloBTagDeepCSV_p71(const Bool_t value) {Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_ = value;}
  void setMu12_DoublePFJets100_CaloBTagDeepCSV_p71(const Bool_t value) {Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_ = value;}
  void setMu12_DoublePFJets200_CaloBTagDeepCSV_p71(const Bool_t value) {Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_ = value;}
  void setMu12_DoublePFJets350_CaloBTagDeepCSV_p71(const Bool_t value) {Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_ = value;}
  void setMu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(const Bool_t value) {Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_ = value;}
  void setMu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(const Bool_t value) {Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_ = value;}
  void setMu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(const Bool_t value) {Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_ = value;}
  void setDoublePFJets40_CaloBTagDeepCSV_p71(const Bool_t value) {DoublePFJets40_CaloBTagDeepCSV_p71_ = value;}
  void setDoublePFJets100_CaloBTagDeepCSV_p71(const Bool_t value) {DoublePFJets100_CaloBTagDeepCSV_p71_ = value;}
  void setDoublePFJets200_CaloBTagDeepCSV_p71(const Bool_t value) {DoublePFJets200_CaloBTagDeepCSV_p71_ = value;}
  void setDoublePFJets350_CaloBTagDeepCSV_p71(const Bool_t value) {DoublePFJets350_CaloBTagDeepCSV_p71_ = value;}
  void setDoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(const Bool_t value) {DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_ = value;}
  void setDoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(const Bool_t value) {DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_ = value;}
  void setPhoton300_NoHE(const Bool_t value) {Photon300_NoHE_ = value;}
  void setMu8_TrkIsoVVL(const Bool_t value) {Mu8_TrkIsoVVL_ = value;}
  void setMu8_DiEle12_CaloIdL_TrackIdL_DZ(const Bool_t value) {Mu8_DiEle12_CaloIdL_TrackIdL_DZ_ = value;}
  void setMu8_DiEle12_CaloIdL_TrackIdL(const Bool_t value) {Mu8_DiEle12_CaloIdL_TrackIdL_ = value;}
  void setMu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ(const Bool_t value) {Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_ = value;}
  void setMu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350(const Bool_t value) {Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_ = value;}
  void setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ(const Bool_t value) {Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_ = value;}
  void setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30(const Bool_t value) {Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_ = value;}
  void setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30(const Bool_t value) {Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_ = value;}
  void setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5(const Bool_t value) {Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_ = value;}
  void setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5(const Bool_t value) {Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_ = value;}
  void setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL(const Bool_t value) {Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_ = value;}
  void setMu17_TrkIsoVVL(const Bool_t value) {Mu17_TrkIsoVVL_ = value;}
  void setMu19_TrkIsoVVL(const Bool_t value) {Mu19_TrkIsoVVL_ = value;}
  void setBTagMu_AK4DiJet20_Mu5(const Bool_t value) {BTagMu_AK4DiJet20_Mu5_ = value;}
  void setBTagMu_AK4DiJet40_Mu5(const Bool_t value) {BTagMu_AK4DiJet40_Mu5_ = value;}
  void setBTagMu_AK4DiJet70_Mu5(const Bool_t value) {BTagMu_AK4DiJet70_Mu5_ = value;}
  void setBTagMu_AK4DiJet110_Mu5(const Bool_t value) {BTagMu_AK4DiJet110_Mu5_ = value;}
  void setBTagMu_AK4DiJet170_Mu5(const Bool_t value) {BTagMu_AK4DiJet170_Mu5_ = value;}
  void setBTagMu_AK4Jet300_Mu5(const Bool_t value) {BTagMu_AK4Jet300_Mu5_ = value;}
  void setBTagMu_AK8DiJet170_Mu5(const Bool_t value) {BTagMu_AK8DiJet170_Mu5_ = value;}
  void setBTagMu_AK8Jet170_DoubleMu5(const Bool_t value) {BTagMu_AK8Jet170_DoubleMu5_ = value;}
  void setBTagMu_AK8Jet300_Mu5(const Bool_t value) {BTagMu_AK8Jet300_Mu5_ = value;}
  void setBTagMu_AK4DiJet20_Mu5_noalgo(const Bool_t value) {BTagMu_AK4DiJet20_Mu5_noalgo_ = value;}
  void setBTagMu_AK4DiJet40_Mu5_noalgo(const Bool_t value) {BTagMu_AK4DiJet40_Mu5_noalgo_ = value;}
  void setBTagMu_AK4DiJet70_Mu5_noalgo(const Bool_t value) {BTagMu_AK4DiJet70_Mu5_noalgo_ = value;}
  void setBTagMu_AK4DiJet110_Mu5_noalgo(const Bool_t value) {BTagMu_AK4DiJet110_Mu5_noalgo_ = value;}
  void setBTagMu_AK4DiJet170_Mu5_noalgo(const Bool_t value) {BTagMu_AK4DiJet170_Mu5_noalgo_ = value;}
  void setBTagMu_AK4Jet300_Mu5_noalgo(const Bool_t value) {BTagMu_AK4Jet300_Mu5_noalgo_ = value;}
  void setBTagMu_AK8DiJet170_Mu5_noalgo(const Bool_t value) {BTagMu_AK8DiJet170_Mu5_noalgo_ = value;}
  void setBTagMu_AK8Jet170_DoubleMu5_noalgo(const Bool_t value) {BTagMu_AK8Jet170_DoubleMu5_noalgo_ = value;}
  void setBTagMu_AK8Jet300_Mu5_noalgo(const Bool_t value) {BTagMu_AK8Jet300_Mu5_noalgo_ = value;}
  void setEle15_Ele8_CaloIdL_TrackIdL_IsoVL(const Bool_t value) {Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_ = value;}
  void setEle23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ(const Bool_t value) {Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_ = value;}
  void setEle23_Ele12_CaloIdL_TrackIdL_IsoVL(const Bool_t value) {Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_ = value;}
  void setMu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ(const Bool_t value) {Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_ = value;}
  void setMu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL(const Bool_t value) {Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_ = value;}
  void setMu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL(const Bool_t value) {Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_ = value;}
  void setMu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ(const Bool_t value) {Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_ = value;}
  void setMu12_DoublePhoton20(const Bool_t value) {Mu12_DoublePhoton20_ = value;}
  void setTriplePhoton_20_20_20_CaloIdLV2(const Bool_t value) {TriplePhoton_20_20_20_CaloIdLV2_ = value;}
  void setTriplePhoton_20_20_20_CaloIdLV2_R9IdVL(const Bool_t value) {TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_ = value;}
  void setTriplePhoton_30_30_10_CaloIdLV2(const Bool_t value) {TriplePhoton_30_30_10_CaloIdLV2_ = value;}
  void setTriplePhoton_30_30_10_CaloIdLV2_R9IdVL(const Bool_t value) {TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_ = value;}
  void setTriplePhoton_35_35_5_CaloIdLV2_R9IdVL(const Bool_t value) {TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_ = value;}
  void setPhoton20(const Bool_t value) {Photon20_ = value;}
  void setPhoton33(const Bool_t value) {Photon33_ = value;}
  void setPhoton50(const Bool_t value) {Photon50_ = value;}
  void setPhoton75(const Bool_t value) {Photon75_ = value;}
  void setPhoton90(const Bool_t value) {Photon90_ = value;}
  void setPhoton120(const Bool_t value) {Photon120_ = value;}
  void setPhoton150(const Bool_t value) {Photon150_ = value;}
  void setPhoton175(const Bool_t value) {Photon175_ = value;}
  void setPhoton200(const Bool_t value) {Photon200_ = value;}
  void setPhoton100EB_TightID_TightIso(const Bool_t value) {Photon100EB_TightID_TightIso_ = value;}
  void setPhoton110EB_TightID_TightIso(const Bool_t value) {Photon110EB_TightID_TightIso_ = value;}
  void setPhoton120EB_TightID_TightIso(const Bool_t value) {Photon120EB_TightID_TightIso_ = value;}
  void setPhoton100EBHE10(const Bool_t value) {Photon100EBHE10_ = value;}
  void setPhoton100EEHE10(const Bool_t value) {Photon100EEHE10_ = value;}
  void setPhoton100EE_TightID_TightIso(const Bool_t value) {Photon100EE_TightID_TightIso_ = value;}
  void setPhoton50_R9Id90_HE10_IsoM(const Bool_t value) {Photon50_R9Id90_HE10_IsoM_ = value;}
  void setPhoton75_R9Id90_HE10_IsoM(const Bool_t value) {Photon75_R9Id90_HE10_IsoM_ = value;}
  void setPhoton75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3(const Bool_t value) {Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_ = value;}
  void setPhoton75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3(const Bool_t value) {Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_ = value;}
  void setPhoton90_R9Id90_HE10_IsoM(const Bool_t value) {Photon90_R9Id90_HE10_IsoM_ = value;}
  void setPhoton120_R9Id90_HE10_IsoM(const Bool_t value) {Photon120_R9Id90_HE10_IsoM_ = value;}
  void setPhoton165_R9Id90_HE10_IsoM(const Bool_t value) {Photon165_R9Id90_HE10_IsoM_ = value;}
  void setPhoton90_CaloIdL_PFHT700(const Bool_t value) {Photon90_CaloIdL_PFHT700_ = value;}
  void setDiphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90(const Bool_t value) {Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_ = value;}
  void setDiphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95(const Bool_t value) {Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_ = value;}
  void setDiphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55(const Bool_t value) {Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_ = value;}
  void setDiphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55(const Bool_t value) {Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_ = value;}
  void setPhoton35_TwoProngs35(const Bool_t value) {Photon35_TwoProngs35_ = value;}
  void setIsoMu24_TwoProngs35(const Bool_t value) {IsoMu24_TwoProngs35_ = value;}
  void setDimuon0_Jpsi_L1_NoOS(const Bool_t value) {Dimuon0_Jpsi_L1_NoOS_ = value;}
  void setDimuon0_Jpsi_NoVertexing_NoOS(const Bool_t value) {Dimuon0_Jpsi_NoVertexing_NoOS_ = value;}
  void setDimuon0_Jpsi(const Bool_t value) {Dimuon0_Jpsi_ = value;}
  void setDimuon0_Jpsi_NoVertexing(const Bool_t value) {Dimuon0_Jpsi_NoVertexing_ = value;}
  void setDimuon0_Jpsi_L1_4R_0er1p5R(const Bool_t value) {Dimuon0_Jpsi_L1_4R_0er1p5R_ = value;}
  void setDimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R(const Bool_t value) {Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_ = value;}
  void setDimuon0_Jpsi3p5_Muon2(const Bool_t value) {Dimuon0_Jpsi3p5_Muon2_ = value;}
  void setDimuon0_Upsilon_L1_4p5(const Bool_t value) {Dimuon0_Upsilon_L1_4p5_ = value;}
  void setDimuon0_Upsilon_L1_5(const Bool_t value) {Dimuon0_Upsilon_L1_5_ = value;}
  void setDimuon0_Upsilon_L1_4p5NoOS(const Bool_t value) {Dimuon0_Upsilon_L1_4p5NoOS_ = value;}
  void setDimuon0_Upsilon_L1_4p5er2p0(const Bool_t value) {Dimuon0_Upsilon_L1_4p5er2p0_ = value;}
  void setDimuon0_Upsilon_L1_4p5er2p0M(const Bool_t value) {Dimuon0_Upsilon_L1_4p5er2p0M_ = value;}
  void setDimuon0_Upsilon_NoVertexing(const Bool_t value) {Dimuon0_Upsilon_NoVertexing_ = value;}
  void setDimuon0_Upsilon_L1_5M(const Bool_t value) {Dimuon0_Upsilon_L1_5M_ = value;}
  void setDimuon0_LowMass_L1_0er1p5R(const Bool_t value) {Dimuon0_LowMass_L1_0er1p5R_ = value;}
  void setDimuon0_LowMass_L1_0er1p5(const Bool_t value) {Dimuon0_LowMass_L1_0er1p5_ = value;}
  void setDimuon0_LowMass(const Bool_t value) {Dimuon0_LowMass_ = value;}
  void setDimuon0_LowMass_L1_4(const Bool_t value) {Dimuon0_LowMass_L1_4_ = value;}
  void setDimuon0_LowMass_L1_4R(const Bool_t value) {Dimuon0_LowMass_L1_4R_ = value;}
  void setDimuon0_LowMass_L1_TM530(const Bool_t value) {Dimuon0_LowMass_L1_TM530_ = value;}
  void setDimuon0_Upsilon_Muon_L1_TM0(const Bool_t value) {Dimuon0_Upsilon_Muon_L1_TM0_ = value;}
  void setDimuon0_Upsilon_Muon_NoL1Mass(const Bool_t value) {Dimuon0_Upsilon_Muon_NoL1Mass_ = value;}
  void setTripleMu_5_3_3_Mass3p8_DZ(const Bool_t value) {TripleMu_5_3_3_Mass3p8_DZ_ = value;}
  void setTripleMu_10_5_5_DZ(const Bool_t value) {TripleMu_10_5_5_DZ_ = value;}
  void setTripleMu_12_10_5(const Bool_t value) {TripleMu_12_10_5_ = value;}
  void setTau3Mu_Mu7_Mu1_TkMu1_Tau15(const Bool_t value) {Tau3Mu_Mu7_Mu1_TkMu1_Tau15_ = value;}
  void setTau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1(const Bool_t value) {Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_ = value;}
  void setTau3Mu_Mu7_Mu1_TkMu1_IsoTau15(const Bool_t value) {Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_ = value;}
  void setTau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1(const Bool_t value) {Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_ = value;}
  void setDoubleMu3_DZ_PFMET50_PFMHT60(const Bool_t value) {DoubleMu3_DZ_PFMET50_PFMHT60_ = value;}
  void setDoubleMu3_DZ_PFMET70_PFMHT70(const Bool_t value) {DoubleMu3_DZ_PFMET70_PFMHT70_ = value;}
  void setDoubleMu3_DZ_PFMET90_PFMHT90(const Bool_t value) {DoubleMu3_DZ_PFMET90_PFMHT90_ = value;}
  void setDoubleMu3_Trk_Tau3mu_NoL1Mass(const Bool_t value) {DoubleMu3_Trk_Tau3mu_NoL1Mass_ = value;}
  void setDoubleMu4_Jpsi_Displaced(const Bool_t value) {DoubleMu4_Jpsi_Displaced_ = value;}
  void setDoubleMu4_Jpsi_NoVertexing(const Bool_t value) {DoubleMu4_Jpsi_NoVertexing_ = value;}
  void setDoubleMu4_JpsiTrkTrk_Displaced(const Bool_t value) {DoubleMu4_JpsiTrkTrk_Displaced_ = value;}
  void setDoubleMu43NoFiltersNoVtx(const Bool_t value) {DoubleMu43NoFiltersNoVtx_ = value;}
  void setDoubleMu48NoFiltersNoVtx(const Bool_t value) {DoubleMu48NoFiltersNoVtx_ = value;}
  void setMu43NoFiltersNoVtx_Photon43_CaloIdL(const Bool_t value) {Mu43NoFiltersNoVtx_Photon43_CaloIdL_ = value;}
  void setMu48NoFiltersNoVtx_Photon48_CaloIdL(const Bool_t value) {Mu48NoFiltersNoVtx_Photon48_CaloIdL_ = value;}
  void setMu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL(const Bool_t value) {Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_ = value;}
  void setMu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL(const Bool_t value) {Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_ = value;}
  void setDoubleMu33NoFiltersNoVtxDisplaced(const Bool_t value) {DoubleMu33NoFiltersNoVtxDisplaced_ = value;}
  void setDoubleMu40NoFiltersNoVtxDisplaced(const Bool_t value) {DoubleMu40NoFiltersNoVtxDisplaced_ = value;}
  void setDoubleMu20_7_Mass0to30_L1_DM4(const Bool_t value) {DoubleMu20_7_Mass0to30_L1_DM4_ = value;}
  void setDoubleMu20_7_Mass0to30_L1_DM4EG(const Bool_t value) {DoubleMu20_7_Mass0to30_L1_DM4EG_ = value;}
  void setHT425(const Bool_t value) {HT425_ = value;}
  void setHT430_DisplacedDijet40_DisplacedTrack(const Bool_t value) {HT430_DisplacedDijet40_DisplacedTrack_ = value;}
  void setHT500_DisplacedDijet40_DisplacedTrack(const Bool_t value) {HT500_DisplacedDijet40_DisplacedTrack_ = value;}
  void setHT430_DisplacedDijet60_DisplacedTrack(const Bool_t value) {HT430_DisplacedDijet60_DisplacedTrack_ = value;}
  void setHT400_DisplacedDijet40_DisplacedTrack(const Bool_t value) {HT400_DisplacedDijet40_DisplacedTrack_ = value;}
  void setHT650_DisplacedDijet60_Inclusive(const Bool_t value) {HT650_DisplacedDijet60_Inclusive_ = value;}
  void setHT550_DisplacedDijet60_Inclusive(const Bool_t value) {HT550_DisplacedDijet60_Inclusive_ = value;}
  void setDiJet110_35_Mjj650_PFMET110(const Bool_t value) {DiJet110_35_Mjj650_PFMET110_ = value;}
  void setDiJet110_35_Mjj650_PFMET120(const Bool_t value) {DiJet110_35_Mjj650_PFMET120_ = value;}
  void setDiJet110_35_Mjj650_PFMET130(const Bool_t value) {DiJet110_35_Mjj650_PFMET130_ = value;}
  void setTripleJet110_35_35_Mjj650_PFMET110(const Bool_t value) {TripleJet110_35_35_Mjj650_PFMET110_ = value;}
  void setTripleJet110_35_35_Mjj650_PFMET120(const Bool_t value) {TripleJet110_35_35_Mjj650_PFMET120_ = value;}
  void setTripleJet110_35_35_Mjj650_PFMET130(const Bool_t value) {TripleJet110_35_35_Mjj650_PFMET130_ = value;}
  void setEle30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned(const Bool_t value) {Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_ = value;}
  void setEle28_eta2p1_WPTight_Gsf_HT150(const Bool_t value) {Ele28_eta2p1_WPTight_Gsf_HT150_ = value;}
  void setEle28_HighEta_SC20_Mass55(const Bool_t value) {Ele28_HighEta_SC20_Mass55_ = value;}
  void setDoubleMu20_7_Mass0to30_Photon23(const Bool_t value) {DoubleMu20_7_Mass0to30_Photon23_ = value;}
  void setEle15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5(const Bool_t value) {Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_ = value;}
  void setEle15_IsoVVVL_PFHT450_PFMET50(const Bool_t value) {Ele15_IsoVVVL_PFHT450_PFMET50_ = value;}
  void setEle15_IsoVVVL_PFHT450(const Bool_t value) {Ele15_IsoVVVL_PFHT450_ = value;}
  void setEle50_IsoVVVL_PFHT450(const Bool_t value) {Ele50_IsoVVVL_PFHT450_ = value;}
  void setEle15_IsoVVVL_PFHT600(const Bool_t value) {Ele15_IsoVVVL_PFHT600_ = value;}
  void setMu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60(const Bool_t value) {Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_ = value;}
  void setMu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60(const Bool_t value) {Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_ = value;}
  void setMu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60(const Bool_t value) {Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_ = value;}
  void setMu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5(const Bool_t value) {Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_ = value;}
  void setMu15_IsoVVVL_PFHT450_PFMET50(const Bool_t value) {Mu15_IsoVVVL_PFHT450_PFMET50_ = value;}
  void setMu15_IsoVVVL_PFHT450(const Bool_t value) {Mu15_IsoVVVL_PFHT450_ = value;}
  void setMu50_IsoVVVL_PFHT450(const Bool_t value) {Mu50_IsoVVVL_PFHT450_ = value;}
  void setMu15_IsoVVVL_PFHT600(const Bool_t value) {Mu15_IsoVVVL_PFHT600_ = value;}
  void setMu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight(const Bool_t value) {Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_ = value;}
  void setMu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight(const Bool_t value) {Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_ = value;}
  void setMu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight(const Bool_t value) {Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_ = value;}
  void setMu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight(const Bool_t value) {Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_ = value;}
  void setMu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight(const Bool_t value) {Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_ = value;}
  void setMu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight(const Bool_t value) {Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_ = value;}
  void setMu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight(const Bool_t value) {Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_ = value;}
  void setMu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight(const Bool_t value) {Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_ = value;}
  void setDimuon10_PsiPrime_Barrel_Seagulls(const Bool_t value) {Dimuon10_PsiPrime_Barrel_Seagulls_ = value;}
  void setDimuon20_Jpsi_Barrel_Seagulls(const Bool_t value) {Dimuon20_Jpsi_Barrel_Seagulls_ = value;}
  void setDimuon12_Upsilon_y1p4(const Bool_t value) {Dimuon12_Upsilon_y1p4_ = value;}
  void setDimuon14_Phi_Barrel_Seagulls(const Bool_t value) {Dimuon14_Phi_Barrel_Seagulls_ = value;}
  void setDimuon18_PsiPrime(const Bool_t value) {Dimuon18_PsiPrime_ = value;}
  void setDimuon25_Jpsi(const Bool_t value) {Dimuon25_Jpsi_ = value;}
  void setDimuon18_PsiPrime_noCorrL1(const Bool_t value) {Dimuon18_PsiPrime_noCorrL1_ = value;}
  void setDimuon24_Upsilon_noCorrL1(const Bool_t value) {Dimuon24_Upsilon_noCorrL1_ = value;}
  void setDimuon24_Phi_noCorrL1(const Bool_t value) {Dimuon24_Phi_noCorrL1_ = value;}
  void setDimuon25_Jpsi_noCorrL1(const Bool_t value) {Dimuon25_Jpsi_noCorrL1_ = value;}
  void setDiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8(const Bool_t value) {DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_ = value;}
  void setDiMu9_Ele9_CaloIdL_TrackIdL_DZ(const Bool_t value) {DiMu9_Ele9_CaloIdL_TrackIdL_DZ_ = value;}
  void setDiMu9_Ele9_CaloIdL_TrackIdL(const Bool_t value) {DiMu9_Ele9_CaloIdL_TrackIdL_ = value;}
  void setDoubleIsoMu20_eta2p1(const Bool_t value) {DoubleIsoMu20_eta2p1_ = value;}
  void setTrkMu12_DoubleTrkMu5NoFiltersNoVtx(const Bool_t value) {TrkMu12_DoubleTrkMu5NoFiltersNoVtx_ = value;}
  void setTrkMu16_DoubleTrkMu6NoFiltersNoVtx(const Bool_t value) {TrkMu16_DoubleTrkMu6NoFiltersNoVtx_ = value;}
  void setTrkMu17_DoubleTrkMu8NoFiltersNoVtx(const Bool_t value) {TrkMu17_DoubleTrkMu8NoFiltersNoVtx_ = value;}
  void setMu8(const Bool_t value) {Mu8_ = value;}
  void setMu17(const Bool_t value) {Mu17_ = value;}
  void setMu19(const Bool_t value) {Mu19_ = value;}
  void setMu17_Photon30_IsoCaloId(const Bool_t value) {Mu17_Photon30_IsoCaloId_ = value;}
  void setEle8_CaloIdL_TrackIdL_IsoVL_PFJet30(const Bool_t value) {Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_ = value;}
  void setEle12_CaloIdL_TrackIdL_IsoVL_PFJet30(const Bool_t value) {Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_ = value;}
  void setEle15_CaloIdL_TrackIdL_IsoVL_PFJet30(const Bool_t value) {Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_ = value;}
  void setEle23_CaloIdL_TrackIdL_IsoVL_PFJet30(const Bool_t value) {Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_ = value;}
  void setEle8_CaloIdM_TrackIdM_PFJet30(const Bool_t value) {Ele8_CaloIdM_TrackIdM_PFJet30_ = value;}
  void setEle17_CaloIdM_TrackIdM_PFJet30(const Bool_t value) {Ele17_CaloIdM_TrackIdM_PFJet30_ = value;}
  void setEle23_CaloIdM_TrackIdM_PFJet30(const Bool_t value) {Ele23_CaloIdM_TrackIdM_PFJet30_ = value;}
  void setEle50_CaloIdVT_GsfTrkIdT_PFJet165(const Bool_t value) {Ele50_CaloIdVT_GsfTrkIdT_PFJet165_ = value;}
  void setEle115_CaloIdVT_GsfTrkIdT(const Bool_t value) {Ele115_CaloIdVT_GsfTrkIdT_ = value;}
  void setEle135_CaloIdVT_GsfTrkIdT(const Bool_t value) {Ele135_CaloIdVT_GsfTrkIdT_ = value;}
  void setEle145_CaloIdVT_GsfTrkIdT(const Bool_t value) {Ele145_CaloIdVT_GsfTrkIdT_ = value;}
  void setEle200_CaloIdVT_GsfTrkIdT(const Bool_t value) {Ele200_CaloIdVT_GsfTrkIdT_ = value;}
  void setEle250_CaloIdVT_GsfTrkIdT(const Bool_t value) {Ele250_CaloIdVT_GsfTrkIdT_ = value;}
  void setEle300_CaloIdVT_GsfTrkIdT(const Bool_t value) {Ele300_CaloIdVT_GsfTrkIdT_ = value;}
  void setPFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5(const Bool_t value) {PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_ = value;}
  void setPFHT330PT30_QuadPFJet_75_60_45_40(const Bool_t value) {PFHT330PT30_QuadPFJet_75_60_45_40_ = value;}
  void setPFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94(const Bool_t value) {PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_ = value;}
  void setPFHT400_SixPFJet32(const Bool_t value) {PFHT400_SixPFJet32_ = value;}
  void setPFHT450_SixPFJet36_PFBTagDeepCSV_1p59(const Bool_t value) {PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_ = value;}
  void setPFHT450_SixPFJet36(const Bool_t value) {PFHT450_SixPFJet36_ = value;}
  void setPFHT350(const Bool_t value) {PFHT350_ = value;}
  void setPFHT350MinPFJet15(const Bool_t value) {PFHT350MinPFJet15_ = value;}
  void setPhoton60_R9Id90_CaloIdL_IsoL(const Bool_t value) {Photon60_R9Id90_CaloIdL_IsoL_ = value;}
  void setPhoton60_R9Id90_CaloIdL_IsoL_DisplacedIdL(const Bool_t value) {Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_ = value;}
  void setPhoton60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15(const Bool_t value) {Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_ = value;}
  void setECALHT800(const Bool_t value) {ECALHT800_ = value;}
  void setDiSC30_18_EIso_AND_HE_Mass70(const Bool_t value) {DiSC30_18_EIso_AND_HE_Mass70_ = value;}
  void setPhysics(const Bool_t value) {Physics_ = value;}
  void setPhysics_part0(const Bool_t value) {Physics_part0_ = value;}
  void setPhysics_part1(const Bool_t value) {Physics_part1_ = value;}
  void setPhysics_part2(const Bool_t value) {Physics_part2_ = value;}
  void setPhysics_part3(const Bool_t value) {Physics_part3_ = value;}
  void setPhysics_part4(const Bool_t value) {Physics_part4_ = value;}
  void setPhysics_part5(const Bool_t value) {Physics_part5_ = value;}
  void setPhysics_part6(const Bool_t value) {Physics_part6_ = value;}
  void setPhysics_part7(const Bool_t value) {Physics_part7_ = value;}
  void setRandom(const Bool_t value) {Random_ = value;}
  void setZeroBias(const Bool_t value) {ZeroBias_ = value;}
  void setZeroBias_Alignment(const Bool_t value) {ZeroBias_Alignment_ = value;}
  void setZeroBias_part0(const Bool_t value) {ZeroBias_part0_ = value;}
  void setZeroBias_part1(const Bool_t value) {ZeroBias_part1_ = value;}
  void setZeroBias_part2(const Bool_t value) {ZeroBias_part2_ = value;}
  void setZeroBias_part3(const Bool_t value) {ZeroBias_part3_ = value;}
  void setZeroBias_part4(const Bool_t value) {ZeroBias_part4_ = value;}
  void setZeroBias_part5(const Bool_t value) {ZeroBias_part5_ = value;}
  void setZeroBias_part6(const Bool_t value) {ZeroBias_part6_ = value;}
  void setZeroBias_part7(const Bool_t value) {ZeroBias_part7_ = value;}
  void setAK4CaloJet30(const Bool_t value) {AK4CaloJet30_ = value;}
  void setAK4CaloJet40(const Bool_t value) {AK4CaloJet40_ = value;}
  void setAK4CaloJet50(const Bool_t value) {AK4CaloJet50_ = value;}
  void setAK4CaloJet80(const Bool_t value) {AK4CaloJet80_ = value;}
  void setAK4CaloJet100(const Bool_t value) {AK4CaloJet100_ = value;}
  void setAK4CaloJet120(const Bool_t value) {AK4CaloJet120_ = value;}
  void setAK4PFJet30(const Bool_t value) {AK4PFJet30_ = value;}
  void setAK4PFJet50(const Bool_t value) {AK4PFJet50_ = value;}
  void setAK4PFJet80(const Bool_t value) {AK4PFJet80_ = value;}
  void setAK4PFJet100(const Bool_t value) {AK4PFJet100_ = value;}
  void setAK4PFJet120(const Bool_t value) {AK4PFJet120_ = value;}
  void setSinglePhoton10_Eta3p1ForPPRef(const Bool_t value) {SinglePhoton10_Eta3p1ForPPRef_ = value;}
  void setSinglePhoton20_Eta3p1ForPPRef(const Bool_t value) {SinglePhoton20_Eta3p1ForPPRef_ = value;}
  void setSinglePhoton30_Eta3p1ForPPRef(const Bool_t value) {SinglePhoton30_Eta3p1ForPPRef_ = value;}
  void setPhoton20_HoverELoose(const Bool_t value) {Photon20_HoverELoose_ = value;}
  void setPhoton30_HoverELoose(const Bool_t value) {Photon30_HoverELoose_ = value;}
  void setEcalCalibration(const Bool_t value) {EcalCalibration_ = value;}
  void setHcalCalibration(const Bool_t value) {HcalCalibration_ = value;}
  void setL1UnpairedBunchBptxMinus(const Bool_t value) {L1UnpairedBunchBptxMinus_ = value;}
  void setL1UnpairedBunchBptxPlus(const Bool_t value) {L1UnpairedBunchBptxPlus_ = value;}
  void setL1NotBptxOR(const Bool_t value) {L1NotBptxOR_ = value;}
  void setL1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142(const Bool_t value) {L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_ = value;}
  void setCDC_L2cosmic_5_er1p0(const Bool_t value) {CDC_L2cosmic_5_er1p0_ = value;}
  void setCDC_L2cosmic_5p5_er1p0(const Bool_t value) {CDC_L2cosmic_5p5_er1p0_ = value;}
  void setHcalNZS(const Bool_t value) {HcalNZS_ = value;}
  void setHcalPhiSym(const Bool_t value) {HcalPhiSym_ = value;}
  void setHcalIsolatedbunch(const Bool_t value) {HcalIsolatedbunch_ = value;}
  void setIsoTrackHB(const Bool_t value) {IsoTrackHB_ = value;}
  void setIsoTrackHE(const Bool_t value) {IsoTrackHE_ = value;}
  void setZeroBias_FirstCollisionAfterAbortGap(const Bool_t value) {ZeroBias_FirstCollisionAfterAbortGap_ = value;}
  void setZeroBias_IsolatedBunches(const Bool_t value) {ZeroBias_IsolatedBunches_ = value;}
  void setZeroBias_FirstCollisionInTrain(const Bool_t value) {ZeroBias_FirstCollisionInTrain_ = value;}
  void setZeroBias_LastCollisionInTrain(const Bool_t value) {ZeroBias_LastCollisionInTrain_ = value;}
  void setZeroBias_FirstBXAfterTrain(const Bool_t value) {ZeroBias_FirstBXAfterTrain_ = value;}
  void setIsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr(const Bool_t value) {IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_ = value;}
  void setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90(const Bool_t value) {MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_ = value;}
  void setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100(const Bool_t value) {MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_ = value;}
  void setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110(const Bool_t value) {MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_ = value;}
  void setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120(const Bool_t value) {MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_ = value;}
  void setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130(const Bool_t value) {MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_ = value;}
  void setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140(const Bool_t value) {MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_ = value;}
  void setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr(const Bool_t value) {MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_ = value;}
  void setMediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr(const Bool_t value) {MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_ = value;}
  void setMediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1(const Bool_t value) {MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_ = value;}
  void setMediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1(const Bool_t value) {MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_ = value;}
  void setMediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1(const Bool_t value) {MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_ = value;}
  void setEle16_Ele12_Ele8_CaloIdL_TrackIdL(const Bool_t value) {Ele16_Ele12_Ele8_CaloIdL_TrackIdL_ = value;}
  void setRsq0p35(const Bool_t value) {Rsq0p35_ = value;}
  void setRsq0p40(const Bool_t value) {Rsq0p40_ = value;}
  void setRsqMR300_Rsq0p09_MR200(const Bool_t value) {RsqMR300_Rsq0p09_MR200_ = value;}
  void setRsqMR320_Rsq0p09_MR200(const Bool_t value) {RsqMR320_Rsq0p09_MR200_ = value;}
  void setRsqMR300_Rsq0p09_MR200_4jet(const Bool_t value) {RsqMR300_Rsq0p09_MR200_4jet_ = value;}
  void setRsqMR320_Rsq0p09_MR200_4jet(const Bool_t value) {RsqMR320_Rsq0p09_MR200_4jet_ = value;}
  void setIsoMu27_MET90(const Bool_t value) {IsoMu27_MET90_ = value;}
  void setDoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg(const Bool_t value) {DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_ = value;}
  void setDoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg(const Bool_t value) {DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_ = value;}
  void setDoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg(const Bool_t value) {DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_ = value;}
  void setDoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg(const Bool_t value) {DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_ = value;}
  void setDoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg(const Bool_t value) {DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_ = value;}
  void setDoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg(const Bool_t value) {DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_ = value;}
  void setDoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg(const Bool_t value) {DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_ = value;}
  void setDoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg(const Bool_t value) {DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_ = value;}
  void setVBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1(const Bool_t value) {VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_ = value;}
  void setVBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1(const Bool_t value) {VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_ = value;}
  void setVBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1(const Bool_t value) {VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_ = value;}
  void setPhoton50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50(const Bool_t value) {Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_ = value;}
  void setPhoton75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3(const Bool_t value) {Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_ = value;}
  void setPhoton75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3(const Bool_t value) {Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_ = value;}
  void setPFMET100_PFMHT100_IDTight_PFHT60(const Bool_t value) {PFMET100_PFMHT100_IDTight_PFHT60_ = value;}
  void setPFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60(const Bool_t value) {PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_ = value;}
  void setPFMETTypeOne100_PFMHT100_IDTight_PFHT60(const Bool_t value) {PFMETTypeOne100_PFMHT100_IDTight_PFHT60_ = value;}
  void setMu18_Mu9_SameSign(const Bool_t value) {Mu18_Mu9_SameSign_ = value;}
  void setMu18_Mu9_SameSign_DZ(const Bool_t value) {Mu18_Mu9_SameSign_DZ_ = value;}
  void setMu18_Mu9(const Bool_t value) {Mu18_Mu9_ = value;}
  void setMu18_Mu9_DZ(const Bool_t value) {Mu18_Mu9_DZ_ = value;}
  void setMu20_Mu10_SameSign(const Bool_t value) {Mu20_Mu10_SameSign_ = value;}
  void setMu20_Mu10_SameSign_DZ(const Bool_t value) {Mu20_Mu10_SameSign_DZ_ = value;}
  void setMu20_Mu10(const Bool_t value) {Mu20_Mu10_ = value;}
  void setMu20_Mu10_DZ(const Bool_t value) {Mu20_Mu10_DZ_ = value;}
  void setMu23_Mu12_SameSign(const Bool_t value) {Mu23_Mu12_SameSign_ = value;}
  void setMu23_Mu12_SameSign_DZ(const Bool_t value) {Mu23_Mu12_SameSign_DZ_ = value;}
  void setMu23_Mu12(const Bool_t value) {Mu23_Mu12_ = value;}
  void setMu23_Mu12_DZ(const Bool_t value) {Mu23_Mu12_DZ_ = value;}
  void setDoubleMu2_Jpsi_DoubleTrk1_Phi1p05(const Bool_t value) {DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_ = value;}
  void setDoubleMu2_Jpsi_DoubleTkMu0_Phi(const Bool_t value) {DoubleMu2_Jpsi_DoubleTkMu0_Phi_ = value;}
  void setDoubleMu3_DCA_PFMET50_PFMHT60(const Bool_t value) {DoubleMu3_DCA_PFMET50_PFMHT60_ = value;}
  void setTripleMu_5_3_3_Mass3p8_DCA(const Bool_t value) {TripleMu_5_3_3_Mass3p8_DCA_ = value;}
  void setQuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1(const Bool_t value) {QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_ = value;}
  void setQuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1(const Bool_t value) {QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_ = value;}
  void setQuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1(const Bool_t value) {QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_ = value;}
  void setQuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2(const Bool_t value) {QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_ = value;}
  void setQuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2(const Bool_t value) {QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_ = value;}
  void setQuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2(const Bool_t value) {QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_ = value;}
  void setQuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2(const Bool_t value) {QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_ = value;}
  void setQuadPFJet98_83_71_15(const Bool_t value) {QuadPFJet98_83_71_15_ = value;}
  void setQuadPFJet103_88_75_15(const Bool_t value) {QuadPFJet103_88_75_15_ = value;}
  void setQuadPFJet105_88_76_15(const Bool_t value) {QuadPFJet105_88_76_15_ = value;}
  void setQuadPFJet111_90_80_15(const Bool_t value) {QuadPFJet111_90_80_15_ = value;}
  void setAK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17(const Bool_t value) {AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_ = value;}
  void setAK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1(const Bool_t value) {AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_ = value;}
  void setAK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02(const Bool_t value) {AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_ = value;}
  void setAK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2(const Bool_t value) {AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_ = value;}
  void setAK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4(const Bool_t value) {AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_ = value;}
  void setDiphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55(const Bool_t value) {Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_ = value;}
  void setDiphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto(const Bool_t value) {Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_ = value;}
  void setMu12_IP6_part0(const Bool_t value) {Mu12_IP6_part0_ = value;}
  void setMu12_IP6_part1(const Bool_t value) {Mu12_IP6_part1_ = value;}
  void setMu12_IP6_part2(const Bool_t value) {Mu12_IP6_part2_ = value;}
  void setMu12_IP6_part3(const Bool_t value) {Mu12_IP6_part3_ = value;}
  void setMu12_IP6_part4(const Bool_t value) {Mu12_IP6_part4_ = value;}
  void setMu9_IP5_part0(const Bool_t value) {Mu9_IP5_part0_ = value;}
  void setMu9_IP5_part1(const Bool_t value) {Mu9_IP5_part1_ = value;}
  void setMu9_IP5_part2(const Bool_t value) {Mu9_IP5_part2_ = value;}
  void setMu9_IP5_part3(const Bool_t value) {Mu9_IP5_part3_ = value;}
  void setMu9_IP5_part4(const Bool_t value) {Mu9_IP5_part4_ = value;}
  void setMu7_IP4_part0(const Bool_t value) {Mu7_IP4_part0_ = value;}
  void setMu7_IP4_part1(const Bool_t value) {Mu7_IP4_part1_ = value;}
  void setMu7_IP4_part2(const Bool_t value) {Mu7_IP4_part2_ = value;}
  void setMu7_IP4_part3(const Bool_t value) {Mu7_IP4_part3_ = value;}
  void setMu7_IP4_part4(const Bool_t value) {Mu7_IP4_part4_ = value;}
  void setMu9_IP4_part0(const Bool_t value) {Mu9_IP4_part0_ = value;}
  void setMu9_IP4_part1(const Bool_t value) {Mu9_IP4_part1_ = value;}
  void setMu9_IP4_part2(const Bool_t value) {Mu9_IP4_part2_ = value;}
  void setMu9_IP4_part3(const Bool_t value) {Mu9_IP4_part3_ = value;}
  void setMu9_IP4_part4(const Bool_t value) {Mu9_IP4_part4_ = value;}
  void setMu8_IP5_part0(const Bool_t value) {Mu8_IP5_part0_ = value;}
  void setMu8_IP5_part1(const Bool_t value) {Mu8_IP5_part1_ = value;}
  void setMu8_IP5_part2(const Bool_t value) {Mu8_IP5_part2_ = value;}
  void setMu8_IP5_part3(const Bool_t value) {Mu8_IP5_part3_ = value;}
  void setMu8_IP5_part4(const Bool_t value) {Mu8_IP5_part4_ = value;}
  void setMu8_IP6_part0(const Bool_t value) {Mu8_IP6_part0_ = value;}
  void setMu8_IP6_part1(const Bool_t value) {Mu8_IP6_part1_ = value;}
  void setMu8_IP6_part2(const Bool_t value) {Mu8_IP6_part2_ = value;}
  void setMu8_IP6_part3(const Bool_t value) {Mu8_IP6_part3_ = value;}
  void setMu8_IP6_part4(const Bool_t value) {Mu8_IP6_part4_ = value;}
  void setMu9_IP6_part0(const Bool_t value) {Mu9_IP6_part0_ = value;}
  void setMu9_IP6_part1(const Bool_t value) {Mu9_IP6_part1_ = value;}
  void setMu9_IP6_part2(const Bool_t value) {Mu9_IP6_part2_ = value;}
  void setMu9_IP6_part3(const Bool_t value) {Mu9_IP6_part3_ = value;}
  void setMu9_IP6_part4(const Bool_t value) {Mu9_IP6_part4_ = value;}
  void setMu8_IP3_part0(const Bool_t value) {Mu8_IP3_part0_ = value;}
  void setMu8_IP3_part1(const Bool_t value) {Mu8_IP3_part1_ = value;}
  void setMu8_IP3_part2(const Bool_t value) {Mu8_IP3_part2_ = value;}
  void setMu8_IP3_part3(const Bool_t value) {Mu8_IP3_part3_ = value;}
  void setMu8_IP3_part4(const Bool_t value) {Mu8_IP3_part4_ = value;}
  void setQuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1(const Bool_t value) {QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_ = value;}
  void setTrkMu6NoFiltersNoVtx(const Bool_t value) {TrkMu6NoFiltersNoVtx_ = value;}
  void setTrkMu16NoFiltersNoVtx(const Bool_t value) {TrkMu16NoFiltersNoVtx_ = value;}
  void setDoubleTrkMu_16_6_NoFiltersNoVtx(const Bool_t value) {DoubleTrkMu_16_6_NoFiltersNoVtx_ = value;}
  void setHLTriggerFinalPath(const Bool_t value) {HLTriggerFinalPath_ = value;}
};

class Svs: public TLorentzVector{
friend class URStreamer;
public:
//  Svs(const Float_t &i_dlen_,const Float_t &i_dlenSig_,const Float_t &i_pAngle_,const Float_t &i_chi2_,const Float_t &i_ndof_,const Float_t &i_x_,const Float_t &i_y_,const Float_t &i_z_):
//    
//  {}
  Svs():
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
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Electrons: public TLorentzVector{
friend class URStreamer;
public:
//  Electrons(const Float_t &i_deltaEtaSC_,const Float_t &i_dr03EcalRecHitSumEt_,const Float_t &i_dr03HcalDepth1TowerSumEt_,const Float_t &i_dr03TkSumPt_,const Float_t &i_dxy_,const Float_t &i_dxyErr_,const Float_t &i_dz_,const Float_t &i_dzErr_,const Float_t &i_eInvMinusPInv_,const Float_t &i_energyErr_,const Float_t &i_hoe_,const Float_t &i_ip3d_,const Float_t &i_miniPFRelIso_all_,const Float_t &i_miniPFRelIso_chg_,const Float_t &i_mvaFall17V1Iso_,const Float_t &i_mvaFall17V1noIso_,const Float_t &i_mvaFall17V2Iso_,const Float_t &i_mvaFall17V2noIso_,const Float_t &i_pfRelIso03_all_,const Float_t &i_pfRelIso03_chg_,const Float_t &i_r9_,const Float_t &i_sieie_,const Float_t &i_sip3d_,const Float_t &i_mvaTTH_,const Int_t &i_charge_,const Int_t &i_cutBased_,const Int_t &i_cutBased_Fall17_V1_,const Int_t &i_jetIdx_,const Int_t &i_pdgId_,const Int_t &i_photonIdx_,const Int_t &i_tightCharge_,const Int_t &i_vidNestedWPBitmap_,const Bool_t &i_convVeto_,const Bool_t &i_cutBased_HEEP_,const Bool_t &i_isPFcand_,const UChar_t &i_lostHits_,const Bool_t &i_mvaFall17V1Iso_WP80_,const Bool_t &i_mvaFall17V1Iso_WP90_,const Bool_t &i_mvaFall17V1Iso_WPL_,const Bool_t &i_mvaFall17V1noIso_WP80_,const Bool_t &i_mvaFall17V1noIso_WP90_,const Bool_t &i_mvaFall17V1noIso_WPL_,const Bool_t &i_mvaFall17V2Iso_WP80_,const Bool_t &i_mvaFall17V2Iso_WP90_,const Bool_t &i_mvaFall17V2Iso_WPL_,const Bool_t &i_mvaFall17V2noIso_WP80_,const Bool_t &i_mvaFall17V2noIso_WP90_,const Bool_t &i_mvaFall17V2noIso_WPL_,const Int_t &i_genPartIdx_,const UChar_t &i_genPartFlav_,const UChar_t &i_cleanmask_):
//    
//  {}
  Electrons():
    TLorentzVector(),
    deltaEtaSC_(0),
    dr03EcalRecHitSumEt_(0),
    dr03HcalDepth1TowerSumEt_(0),
    dr03TkSumPt_(0),
    dxy_(0),
    dxyErr_(0),
    dz_(0),
    dzErr_(0),
    eInvMinusPInv_(0),
    energyErr_(0),
    hoe_(0),
    ip3d_(0),
    miniPFRelIso_all_(0),
    miniPFRelIso_chg_(0),
    mvaFall17V1Iso_(0),
    mvaFall17V1noIso_(0),
    mvaFall17V2Iso_(0),
    mvaFall17V2noIso_(0),
    pfRelIso03_all_(0),
    pfRelIso03_chg_(0),
    r9_(0),
    sieie_(0),
    sip3d_(0),
    mvaTTH_(0),
    charge_(0),
    cutBased_(0),
    cutBased_Fall17_V1_(0),
    jetIdx_(0),
    pdgId_(0),
    photonIdx_(0),
    tightCharge_(0),
    vidNestedWPBitmap_(0),
    convVeto_(0),
    cutBased_HEEP_(0),
    isPFcand_(0),
    lostHits_(0),
    mvaFall17V1Iso_WP80_(0),
    mvaFall17V1Iso_WP90_(0),
    mvaFall17V1Iso_WPL_(0),
    mvaFall17V1noIso_WP80_(0),
    mvaFall17V1noIso_WP90_(0),
    mvaFall17V1noIso_WPL_(0),
    mvaFall17V2Iso_WP80_(0),
    mvaFall17V2Iso_WP90_(0),
    mvaFall17V2Iso_WPL_(0),
    mvaFall17V2noIso_WP80_(0),
    mvaFall17V2noIso_WP90_(0),
    mvaFall17V2noIso_WPL_(0),
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
  Float_t eInvMinusPInv() const {return eInvMinusPInv_;}
  Float_t energyErr() const {return energyErr_;}
  Float_t hoe() const {return hoe_;}
  Float_t ip3d() const {return ip3d_;}
  Float_t miniPFRelIso_all() const {return miniPFRelIso_all_;}
  Float_t miniPFRelIso_chg() const {return miniPFRelIso_chg_;}
  Float_t mvaFall17V1Iso() const {return mvaFall17V1Iso_;}
  Float_t mvaFall17V1noIso() const {return mvaFall17V1noIso_;}
  Float_t mvaFall17V2Iso() const {return mvaFall17V2Iso_;}
  Float_t mvaFall17V2noIso() const {return mvaFall17V2noIso_;}
  Float_t pfRelIso03_all() const {return pfRelIso03_all_;}
  Float_t pfRelIso03_chg() const {return pfRelIso03_chg_;}
  Float_t r9() const {return r9_;}
  Float_t sieie() const {return sieie_;}
  Float_t sip3d() const {return sip3d_;}
  Float_t mvaTTH() const {return mvaTTH_;}
  Int_t charge() const {return charge_;}
  Int_t cutBased() const {return cutBased_;}
  Int_t cutBased_Fall17_V1() const {return cutBased_Fall17_V1_;}
  Int_t jetIdx() const {return jetIdx_;}
  Int_t pdgId() const {return pdgId_;}
  Int_t photonIdx() const {return photonIdx_;}
  Int_t tightCharge() const {return tightCharge_;}
  Int_t vidNestedWPBitmap() const {return vidNestedWPBitmap_;}
  Bool_t convVeto() const {return convVeto_;}
  Bool_t cutBased_HEEP() const {return cutBased_HEEP_;}
  Bool_t isPFcand() const {return isPFcand_;}
  UChar_t lostHits() const {return lostHits_;}
  Bool_t mvaFall17V1Iso_WP80() const {return mvaFall17V1Iso_WP80_;}
  Bool_t mvaFall17V1Iso_WP90() const {return mvaFall17V1Iso_WP90_;}
  Bool_t mvaFall17V1Iso_WPL() const {return mvaFall17V1Iso_WPL_;}
  Bool_t mvaFall17V1noIso_WP80() const {return mvaFall17V1noIso_WP80_;}
  Bool_t mvaFall17V1noIso_WP90() const {return mvaFall17V1noIso_WP90_;}
  Bool_t mvaFall17V1noIso_WPL() const {return mvaFall17V1noIso_WPL_;}
  Bool_t mvaFall17V2Iso_WP80() const {return mvaFall17V2Iso_WP80_;}
  Bool_t mvaFall17V2Iso_WP90() const {return mvaFall17V2Iso_WP90_;}
  Bool_t mvaFall17V2Iso_WPL() const {return mvaFall17V2Iso_WPL_;}
  Bool_t mvaFall17V2noIso_WP80() const {return mvaFall17V2noIso_WP80_;}
  Bool_t mvaFall17V2noIso_WP90() const {return mvaFall17V2noIso_WP90_;}
  Bool_t mvaFall17V2noIso_WPL() const {return mvaFall17V2noIso_WPL_;}
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
  Float_t eInvMinusPInv_;
  Float_t energyErr_;
  Float_t hoe_;
  Float_t ip3d_;
  Float_t miniPFRelIso_all_;
  Float_t miniPFRelIso_chg_;
  Float_t mvaFall17V1Iso_;
  Float_t mvaFall17V1noIso_;
  Float_t mvaFall17V2Iso_;
  Float_t mvaFall17V2noIso_;
  Float_t pfRelIso03_all_;
  Float_t pfRelIso03_chg_;
  Float_t r9_;
  Float_t sieie_;
  Float_t sip3d_;
  Float_t mvaTTH_;
  Int_t charge_;
  Int_t cutBased_;
  Int_t cutBased_Fall17_V1_;
  Int_t jetIdx_;
  Int_t pdgId_;
  Int_t photonIdx_;
  Int_t tightCharge_;
  Int_t vidNestedWPBitmap_;
  Bool_t convVeto_;
  Bool_t cutBased_HEEP_;
  Bool_t isPFcand_;
  UChar_t lostHits_;
  Bool_t mvaFall17V1Iso_WP80_;
  Bool_t mvaFall17V1Iso_WP90_;
  Bool_t mvaFall17V1Iso_WPL_;
  Bool_t mvaFall17V1noIso_WP80_;
  Bool_t mvaFall17V1noIso_WP90_;
  Bool_t mvaFall17V1noIso_WPL_;
  Bool_t mvaFall17V2Iso_WP80_;
  Bool_t mvaFall17V2Iso_WP90_;
  Bool_t mvaFall17V2Iso_WPL_;
  Bool_t mvaFall17V2noIso_WP80_;
  Bool_t mvaFall17V2noIso_WP90_;
  Bool_t mvaFall17V2noIso_WPL_;
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
  void seteInvMinusPInv(const Float_t value) {eInvMinusPInv_ = value;}
  void setenergyErr(const Float_t value) {energyErr_ = value;}
  void sethoe(const Float_t value) {hoe_ = value;}
  void setip3d(const Float_t value) {ip3d_ = value;}
  void setminiPFRelIso_all(const Float_t value) {miniPFRelIso_all_ = value;}
  void setminiPFRelIso_chg(const Float_t value) {miniPFRelIso_chg_ = value;}
  void setmvaFall17V1Iso(const Float_t value) {mvaFall17V1Iso_ = value;}
  void setmvaFall17V1noIso(const Float_t value) {mvaFall17V1noIso_ = value;}
  void setmvaFall17V2Iso(const Float_t value) {mvaFall17V2Iso_ = value;}
  void setmvaFall17V2noIso(const Float_t value) {mvaFall17V2noIso_ = value;}
  void setpfRelIso03_all(const Float_t value) {pfRelIso03_all_ = value;}
  void setpfRelIso03_chg(const Float_t value) {pfRelIso03_chg_ = value;}
  void setr9(const Float_t value) {r9_ = value;}
  void setsieie(const Float_t value) {sieie_ = value;}
  void setsip3d(const Float_t value) {sip3d_ = value;}
  void setmvaTTH(const Float_t value) {mvaTTH_ = value;}
  void setcharge(const Int_t value) {charge_ = value;}
  void setcutBased(const Int_t value) {cutBased_ = value;}
  void setcutBased_Fall17_V1(const Int_t value) {cutBased_Fall17_V1_ = value;}
  void setjetIdx(const Int_t value) {jetIdx_ = value;}
  void setpdgId(const Int_t value) {pdgId_ = value;}
  void setphotonIdx(const Int_t value) {photonIdx_ = value;}
  void settightCharge(const Int_t value) {tightCharge_ = value;}
  void setvidNestedWPBitmap(const Int_t value) {vidNestedWPBitmap_ = value;}
  void setconvVeto(const Bool_t value) {convVeto_ = value;}
  void setcutBased_HEEP(const Bool_t value) {cutBased_HEEP_ = value;}
  void setisPFcand(const Bool_t value) {isPFcand_ = value;}
  void setlostHits(const UChar_t value) {lostHits_ = value;}
  void setmvaFall17V1Iso_WP80(const Bool_t value) {mvaFall17V1Iso_WP80_ = value;}
  void setmvaFall17V1Iso_WP90(const Bool_t value) {mvaFall17V1Iso_WP90_ = value;}
  void setmvaFall17V1Iso_WPL(const Bool_t value) {mvaFall17V1Iso_WPL_ = value;}
  void setmvaFall17V1noIso_WP80(const Bool_t value) {mvaFall17V1noIso_WP80_ = value;}
  void setmvaFall17V1noIso_WP90(const Bool_t value) {mvaFall17V1noIso_WP90_ = value;}
  void setmvaFall17V1noIso_WPL(const Bool_t value) {mvaFall17V1noIso_WPL_ = value;}
  void setmvaFall17V2Iso_WP80(const Bool_t value) {mvaFall17V2Iso_WP80_ = value;}
  void setmvaFall17V2Iso_WP90(const Bool_t value) {mvaFall17V2Iso_WP90_ = value;}
  void setmvaFall17V2Iso_WPL(const Bool_t value) {mvaFall17V2Iso_WPL_ = value;}
  void setmvaFall17V2noIso_WP80(const Bool_t value) {mvaFall17V2noIso_WP80_ = value;}
  void setmvaFall17V2noIso_WP90(const Bool_t value) {mvaFall17V2noIso_WP90_ = value;}
  void setmvaFall17V2noIso_WPL(const Bool_t value) {mvaFall17V2noIso_WPL_ = value;}
  void setgenPartIdx(const Int_t value) {genPartIdx_ = value;}
  void setgenPartFlav(const UChar_t value) {genPartFlav_ = value;}
  void setcleanmask(const UChar_t value) {cleanmask_ = value;}
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
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
  void setLorentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Genmet: public TLorentzVector{
friend class URStreamer;
public:
//  Genmet():
//    
//  {}
  Genmet():
    TLorentzVector()
    
  {}
  
private:
  
  
  void setLorentzVector(float pt, float phi){SetPtEtaPhiM(pt, 0., phi, 0.);}
};

class Fatjets: public TLorentzVector{
friend class URStreamer;
public:
//  Fatjets(const Float_t &i_area_,const Float_t &i_btagCMVA_,const Float_t &i_btagCSVV2_,const Float_t &i_btagDeepB_,const Float_t &i_btagHbb_,const Float_t &i_msoftdrop_,const Float_t &i_n2b1_,const Float_t &i_n3b1_,const Float_t &i_rawFactor_,const Float_t &i_tau1_,const Float_t &i_tau2_,const Float_t &i_tau3_,const Float_t &i_tau4_,const Int_t &i_jetId_,const Int_t &i_subJetIdx1_,const Int_t &i_subJetIdx2_):
//    
//  {}
  Fatjets():
    TLorentzVector(),
    area_(0),
    btagCMVA_(0),
    btagCSVV2_(0),
    btagDeepB_(0),
    btagHbb_(0),
    msoftdrop_(0),
    n2b1_(0),
    n3b1_(0),
    rawFactor_(0),
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
  Float_t n2b1() const {return n2b1_;}
  Float_t n3b1() const {return n3b1_;}
  Float_t rawFactor() const {return rawFactor_;}
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
  Float_t n2b1_;
  Float_t n3b1_;
  Float_t rawFactor_;
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
  void setn2b1(const Float_t value) {n2b1_ = value;}
  void setn3b1(const Float_t value) {n3b1_ = value;}
  void setrawFactor(const Float_t value) {rawFactor_ = value;}
  void settau1(const Float_t value) {tau1_ = value;}
  void settau2(const Float_t value) {tau2_ = value;}
  void settau3(const Float_t value) {tau3_ = value;}
  void settau4(const Float_t value) {tau4_ = value;}
  void setjetId(const Int_t value) {jetId_ = value;}
  void setsubJetIdx1(const Int_t value) {subJetIdx1_ = value;}
  void setsubJetIdx2(const Int_t value) {subJetIdx2_ = value;}
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Subjets: public TLorentzVector{
friend class URStreamer;
public:
//  Subjets(const Float_t &i_btagCMVA_,const Float_t &i_btagCSVV2_,const Float_t &i_btagDeepB_,const Float_t &i_n2b1_,const Float_t &i_n3b1_,const Float_t &i_tau1_,const Float_t &i_tau2_,const Float_t &i_tau3_,const Float_t &i_tau4_):
//    
//  {}
  Subjets():
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
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
};

class Subgenjetak8s: public TLorentzVector{
friend class URStreamer;
public:
//  Subgenjetak8s():
//    
//  {}
  Subgenjetak8s():
    TLorentzVector()
    
  {}
  
private:
  
  
  void setLorentzVector(float pt, float eta, float phi, float mass){SetPtEtaPhiM(pt, eta, phi, mass);}
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
  UInt_t nSubGenJetAK8;
  UInt_t nGenVisTau;
  Float_t genWeight;
  UInt_t nLHEPdfWeight;
  vector<Float_t> LHEPdfWeight;
  UInt_t nLHEScaleWeight;
  vector<Float_t> LHEScaleWeight;
  UInt_t nPSWeight;
  vector<Float_t> PSWeight;
  UInt_t nIsoTrack;
  UInt_t nJet;
  UInt_t nLHEPart;
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
    btagWeight_CSVV2_(0),
    btagWeight_DeepCSVB_(0),
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
    Electron_eInvMinusPInv_(0),
    Electron_energyErr_(0),
    Electron_eta_(0),
    Electron_hoe_(0),
    Electron_ip3d_(0),
    Electron_mass_(0),
    Electron_miniPFRelIso_all_(0),
    Electron_miniPFRelIso_chg_(0),
    Electron_mvaFall17V1Iso_(0),
    Electron_mvaFall17V1noIso_(0),
    Electron_mvaFall17V2Iso_(0),
    Electron_mvaFall17V2noIso_(0),
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
    Electron_cutBased_Fall17_V1_(0),
    Electron_jetIdx_(0),
    Electron_pdgId_(0),
    Electron_photonIdx_(0),
    Electron_tightCharge_(0),
    Electron_vidNestedWPBitmap_(0),
    Electron_convVeto_(0),
    Electron_cutBased_HEEP_(0),
    Electron_isPFcand_(0),
    Electron_lostHits_(0),
    Electron_mvaFall17V1Iso_WP80_(0),
    Electron_mvaFall17V1Iso_WP90_(0),
    Electron_mvaFall17V1Iso_WPL_(0),
    Electron_mvaFall17V1noIso_WP80_(0),
    Electron_mvaFall17V1noIso_WP90_(0),
    Electron_mvaFall17V1noIso_WPL_(0),
    Electron_mvaFall17V2Iso_WP80_(0),
    Electron_mvaFall17V2Iso_WP90_(0),
    Electron_mvaFall17V2Iso_WPL_(0),
    Electron_mvaFall17V2noIso_WP80_(0),
    Electron_mvaFall17V2noIso_WP90_(0),
    Electron_mvaFall17V2noIso_WPL_(0),
    nFatJet(0),
    FatJet_area_(0),
    FatJet_btagCMVA_(0),
    FatJet_btagCSVV2_(0),
    FatJet_btagDeepB_(0),
    FatJet_btagHbb_(0),
    FatJet_eta_(0),
    FatJet_mass_(0),
    FatJet_msoftdrop_(0),
    FatJet_n2b1_(0),
    FatJet_n3b1_(0),
    FatJet_phi_(0),
    FatJet_pt_(0),
    FatJet_rawFactor_(0),
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
    nSubGenJetAK8(0),
    SubGenJetAK8_eta_(0),
    SubGenJetAK8_mass_(0),
    SubGenJetAK8_phi_(0),
    SubGenJetAK8_pt_(0),
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
    nPSWeight(0),
    PSWeight(0),
    nIsoTrack(0),
    IsoTrack_dxy_(0),
    IsoTrack_dz_(0),
    IsoTrack_eta_(0),
    IsoTrack_pfRelIso03_all_(0),
    IsoTrack_pfRelIso03_chg_(0),
    IsoTrack_phi_(0),
    IsoTrack_pt_(0),
    IsoTrack_miniPFRelIso_all_(0),
    IsoTrack_miniPFRelIso_chg_(0),
    IsoTrack_pdgId_(0),
    IsoTrack_isHighPurityTrack_(0),
    IsoTrack_isPFcand_(0),
    nJet(0),
    Jet_area_(0),
    Jet_btagCMVA_(0),
    Jet_btagCSVV2_(0),
    Jet_btagDeepB_(0),
    Jet_btagDeepC_(0),
    Jet_btagDeepFlavB_(0),
    Jet_chEmEF_(0),
    Jet_chHEF_(0),
    Jet_eta_(0),
    Jet_mass_(0),
    Jet_muEF_(0),
    Jet_neEmEF_(0),
    Jet_neHEF_(0),
    Jet_phi_(0),
    Jet_pt_(0),
    Jet_qgl_(0),
    Jet_rawFactor_(0),
    Jet_bRegCorr_(0),
    Jet_bRegRes_(0),
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
    nLHEPart(0),
    LHEPart_pt_(0),
    LHEPart_eta_(0),
    LHEPart_phi_(0),
    LHEPart_mass_(0),
    LHEPart_pdgId_(0),
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
    Photon_cutBasedBitmap_(0),
    Photon_electronIdx_(0),
    Photon_jetIdx_(0),
    Photon_pdgId_(0),
    Photon_vidNestedWPBitmap_(0),
    Photon_electronVeto_(0),
    Photon_isScEtaEB_(0),
    Photon_isScEtaEE_(0),
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
    Tau_rawIsodR03_(0),
    Tau_charge_(0),
    Tau_decayMode_(0),
    Tau_jetIdx_(0),
    Tau_rawAntiEleCat_(0),
    Tau_idAntiEle_(0),
    Tau_idAntiMu_(0),
    Tau_idDecayMode_(0),
    Tau_idDecayModeNewDMs_(0),
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
    HLT_AK8PFJet380_TrimMass30_(0),
    HLT_AK8PFJet400_TrimMass30_(0),
    HLT_AK8PFJet420_TrimMass30_(0),
    HLT_AK8PFHT750_TrimMass50_(0),
    HLT_AK8PFHT800_TrimMass50_(0),
    HLT_AK8PFHT850_TrimMass50_(0),
    HLT_AK8PFHT900_TrimMass50_(0),
    HLT_CaloJet500_NoJetID_(0),
    HLT_CaloJet550_NoJetID_(0),
    HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_(0),
    HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_(0),
    HLT_Trimuon5_3p5_2_Upsilon_Muon_(0),
    HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_(0),
    HLT_DoubleEle25_CaloIdL_MW_(0),
    HLT_DoubleEle27_CaloIdL_MW_(0),
    HLT_DoubleEle33_CaloIdL_MW_(0),
    HLT_DoubleEle24_eta2p1_WPTight_Gsf_(0),
    HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_(0),
    HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_(0),
    HLT_Ele27_Ele37_CaloIdL_MW_(0),
    HLT_Mu27_Ele37_CaloIdL_MW_(0),
    HLT_Mu37_Ele27_CaloIdL_MW_(0),
    HLT_Mu37_TkMu27_(0),
    HLT_DoubleMu4_3_Bs_(0),
    HLT_DoubleMu4_3_Jpsi_(0),
    HLT_DoubleMu4_JpsiTrk_Displaced_(0),
    HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_(0),
    HLT_DoubleMu3_Trk_Tau3mu_(0),
    HLT_DoubleMu3_TkMu_DsTau3Mu_(0),
    HLT_DoubleMu4_PsiPrimeTrk_Displaced_(0),
    HLT_DoubleMu4_Mass3p8_DZ_PFHT350_(0),
    HLT_Mu3_PFJet40_(0),
    HLT_Mu7p5_L2Mu2_Jpsi_(0),
    HLT_Mu7p5_L2Mu2_Upsilon_(0),
    HLT_Mu7p5_Track2_Jpsi_(0),
    HLT_Mu7p5_Track3p5_Jpsi_(0),
    HLT_Mu7p5_Track7_Jpsi_(0),
    HLT_Mu7p5_Track2_Upsilon_(0),
    HLT_Mu7p5_Track3p5_Upsilon_(0),
    HLT_Mu7p5_Track7_Upsilon_(0),
    HLT_Mu3_L1SingleMu5orSingleMu7_(0),
    HLT_DoublePhoton33_CaloIdL_(0),
    HLT_DoublePhoton70_(0),
    HLT_DoublePhoton85_(0),
    HLT_Ele20_WPTight_Gsf_(0),
    HLT_Ele15_WPLoose_Gsf_(0),
    HLT_Ele17_WPLoose_Gsf_(0),
    HLT_Ele20_WPLoose_Gsf_(0),
    HLT_Ele20_eta2p1_WPLoose_Gsf_(0),
    HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_(0),
    HLT_Ele27_WPTight_Gsf_(0),
    HLT_Ele28_WPTight_Gsf_(0),
    HLT_Ele30_WPTight_Gsf_(0),
    HLT_Ele32_WPTight_Gsf_(0),
    HLT_Ele35_WPTight_Gsf_(0),
    HLT_Ele35_WPTight_Gsf_L1EGMT_(0),
    HLT_Ele38_WPTight_Gsf_(0),
    HLT_Ele40_WPTight_Gsf_(0),
    HLT_Ele32_WPTight_Gsf_L1DoubleEG_(0),
    HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_(0),
    HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_(0),
    HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_(0),
    HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_(0),
    HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_(0),
    HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_(0),
    HLT_HT450_Beamspot_(0),
    HLT_HT300_Beamspot_(0),
    HLT_ZeroBias_Beamspot_(0),
    HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_(0),
    HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_(0),
    HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_(0),
    HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_(0),
    HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_(0),
    HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_(0),
    HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_(0),
    HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_(0),
    HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_(0),
    HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_(0),
    HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_(0),
    HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_(0),
    HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_(0),
    HLT_IsoMu20_(0),
    HLT_IsoMu24_(0),
    HLT_IsoMu24_eta2p1_(0),
    HLT_IsoMu27_(0),
    HLT_IsoMu30_(0),
    HLT_UncorrectedJetE30_NoBPTX_(0),
    HLT_UncorrectedJetE30_NoBPTX3BX_(0),
    HLT_UncorrectedJetE60_NoBPTX3BX_(0),
    HLT_UncorrectedJetE70_NoBPTX3BX_(0),
    HLT_L1SingleMu18_(0),
    HLT_L1SingleMu25_(0),
    HLT_L2Mu10_(0),
    HLT_L2Mu10_NoVertex_NoBPTX3BX_(0),
    HLT_L2Mu10_NoVertex_NoBPTX_(0),
    HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_(0),
    HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_(0),
    HLT_L2Mu50_(0),
    HLT_L2Mu23NoVtx_2Cha_(0),
    HLT_L2Mu23NoVtx_2Cha_CosmicSeed_(0),
    HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_(0),
    HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_(0),
    HLT_DoubleL2Mu50_(0),
    HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_(0),
    HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_(0),
    HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_(0),
    HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_(0),
    HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_(0),
    HLT_DoubleL2Mu23NoVtx_2Cha_(0),
    HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_(0),
    HLT_DoubleL2Mu25NoVtx_2Cha_(0),
    HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_(0),
    HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_(0),
    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_(0),
    HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_(0),
    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_(0),
    HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_(0),
    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_(0),
    HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_(0),
    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_(0),
    HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_(0),
    HLT_Mu25_TkMu0_Onia_(0),
    HLT_Mu30_TkMu0_Psi_(0),
    HLT_Mu30_TkMu0_Upsilon_(0),
    HLT_Mu20_TkMu0_Phi_(0),
    HLT_Mu25_TkMu0_Phi_(0),
    HLT_Mu12_(0),
    HLT_Mu15_(0),
    HLT_Mu20_(0),
    HLT_Mu27_(0),
    HLT_Mu50_(0),
    HLT_Mu55_(0),
    HLT_OldMu100_(0),
    HLT_TkMu100_(0),
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
    HLT_AK8PFJet15_(0),
    HLT_AK8PFJet25_(0),
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
    HLT_AK8PFJet550_(0),
    HLT_PFJet15_(0),
    HLT_PFJet25_(0),
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
    HLT_PFJet550_(0),
    HLT_PFJetFwd15_(0),
    HLT_PFJetFwd25_(0),
    HLT_PFJetFwd40_(0),
    HLT_PFJetFwd60_(0),
    HLT_PFJetFwd80_(0),
    HLT_PFJetFwd140_(0),
    HLT_PFJetFwd200_(0),
    HLT_PFJetFwd260_(0),
    HLT_PFJetFwd320_(0),
    HLT_PFJetFwd400_(0),
    HLT_PFJetFwd450_(0),
    HLT_PFJetFwd500_(0),
    HLT_AK8PFJetFwd15_(0),
    HLT_AK8PFJetFwd25_(0),
    HLT_AK8PFJetFwd40_(0),
    HLT_AK8PFJetFwd60_(0),
    HLT_AK8PFJetFwd80_(0),
    HLT_AK8PFJetFwd140_(0),
    HLT_AK8PFJetFwd200_(0),
    HLT_AK8PFJetFwd260_(0),
    HLT_AK8PFJetFwd320_(0),
    HLT_AK8PFJetFwd400_(0),
    HLT_AK8PFJetFwd450_(0),
    HLT_AK8PFJetFwd500_(0),
    HLT_PFHT180_(0),
    HLT_PFHT250_(0),
    HLT_PFHT370_(0),
    HLT_PFHT430_(0),
    HLT_PFHT510_(0),
    HLT_PFHT590_(0),
    HLT_PFHT680_(0),
    HLT_PFHT780_(0),
    HLT_PFHT890_(0),
    HLT_PFHT1050_(0),
    HLT_PFHT500_PFMET100_PFMHT100_IDTight_(0),
    HLT_PFHT500_PFMET110_PFMHT110_IDTight_(0),
    HLT_PFHT700_PFMET85_PFMHT85_IDTight_(0),
    HLT_PFHT700_PFMET95_PFMHT95_IDTight_(0),
    HLT_PFHT800_PFMET75_PFMHT75_IDTight_(0),
    HLT_PFHT800_PFMET85_PFMHT85_IDTight_(0),
    HLT_PFMET110_PFMHT110_IDTight_(0),
    HLT_PFMET120_PFMHT120_IDTight_(0),
    HLT_PFMET130_PFMHT130_IDTight_(0),
    HLT_PFMET140_PFMHT140_IDTight_(0),
    HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_(0),
    HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_(0),
    HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_(0),
    HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_(0),
    HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_(0),
    HLT_PFMET120_PFMHT120_IDTight_PFHT60_(0),
    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_(0),
    HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_(0),
    HLT_PFMETTypeOne110_PFMHT110_IDTight_(0),
    HLT_PFMETTypeOne120_PFMHT120_IDTight_(0),
    HLT_PFMETTypeOne130_PFMHT130_IDTight_(0),
    HLT_PFMETTypeOne140_PFMHT140_IDTight_(0),
    HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_(0),
    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_(0),
    HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_(0),
    HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_(0),
    HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_(0),
    HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_(0),
    HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_(0),
    HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_(0),
    HLT_L1ETMHadSeeds_(0),
    HLT_CaloMHT90_(0),
    HLT_CaloMET80_NotCleaned_(0),
    HLT_CaloMET90_NotCleaned_(0),
    HLT_CaloMET100_NotCleaned_(0),
    HLT_CaloMET110_NotCleaned_(0),
    HLT_CaloMET250_NotCleaned_(0),
    HLT_CaloMET70_HBHECleaned_(0),
    HLT_CaloMET80_HBHECleaned_(0),
    HLT_CaloMET90_HBHECleaned_(0),
    HLT_CaloMET100_HBHECleaned_(0),
    HLT_CaloMET250_HBHECleaned_(0),
    HLT_CaloMET300_HBHECleaned_(0),
    HLT_CaloMET350_HBHECleaned_(0),
    HLT_PFMET200_NotCleaned_(0),
    HLT_PFMET200_HBHECleaned_(0),
    HLT_PFMET250_HBHECleaned_(0),
    HLT_PFMET300_HBHECleaned_(0),
    HLT_PFMET200_HBHE_BeamHaloCleaned_(0),
    HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_(0),
    HLT_MET105_IsoTrk50_(0),
    HLT_MET120_IsoTrk50_(0),
    HLT_SingleJet30_Mu12_SinglePFJet40_(0),
    HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_(0),
    HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_(0),
    HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_(0),
    HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_(0),
    HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    HLT_DoublePFJets40_CaloBTagDeepCSV_p71_(0),
    HLT_DoublePFJets100_CaloBTagDeepCSV_p71_(0),
    HLT_DoublePFJets200_CaloBTagDeepCSV_p71_(0),
    HLT_DoublePFJets350_CaloBTagDeepCSV_p71_(0),
    HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_(0),
    HLT_Photon300_NoHE_(0),
    HLT_Mu8_TrkIsoVVL_(0),
    HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_(0),
    HLT_Mu8_DiEle12_CaloIdL_TrackIdL_(0),
    HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_(0),
    HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_(0),
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_(0),
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_(0),
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_(0),
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_(0),
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu17_TrkIsoVVL_(0),
    HLT_Mu19_TrkIsoVVL_(0),
    HLT_BTagMu_AK4DiJet20_Mu5_(0),
    HLT_BTagMu_AK4DiJet40_Mu5_(0),
    HLT_BTagMu_AK4DiJet70_Mu5_(0),
    HLT_BTagMu_AK4DiJet110_Mu5_(0),
    HLT_BTagMu_AK4DiJet170_Mu5_(0),
    HLT_BTagMu_AK4Jet300_Mu5_(0),
    HLT_BTagMu_AK8DiJet170_Mu5_(0),
    HLT_BTagMu_AK8Jet170_DoubleMu5_(0),
    HLT_BTagMu_AK8Jet300_Mu5_(0),
    HLT_BTagMu_AK4DiJet20_Mu5_noalgo_(0),
    HLT_BTagMu_AK4DiJet40_Mu5_noalgo_(0),
    HLT_BTagMu_AK4DiJet70_Mu5_noalgo_(0),
    HLT_BTagMu_AK4DiJet110_Mu5_noalgo_(0),
    HLT_BTagMu_AK4DiJet170_Mu5_noalgo_(0),
    HLT_BTagMu_AK4Jet300_Mu5_noalgo_(0),
    HLT_BTagMu_AK8DiJet170_Mu5_noalgo_(0),
    HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_(0),
    HLT_BTagMu_AK8Jet300_Mu5_noalgo_(0),
    HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(0),
    HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_(0),
    HLT_Mu12_DoublePhoton20_(0),
    HLT_TriplePhoton_20_20_20_CaloIdLV2_(0),
    HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_(0),
    HLT_TriplePhoton_30_30_10_CaloIdLV2_(0),
    HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_(0),
    HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_(0),
    HLT_Photon20_(0),
    HLT_Photon33_(0),
    HLT_Photon50_(0),
    HLT_Photon75_(0),
    HLT_Photon90_(0),
    HLT_Photon120_(0),
    HLT_Photon150_(0),
    HLT_Photon175_(0),
    HLT_Photon200_(0),
    HLT_Photon100EB_TightID_TightIso_(0),
    HLT_Photon110EB_TightID_TightIso_(0),
    HLT_Photon120EB_TightID_TightIso_(0),
    HLT_Photon100EBHE10_(0),
    HLT_Photon100EEHE10_(0),
    HLT_Photon100EE_TightID_TightIso_(0),
    HLT_Photon50_R9Id90_HE10_IsoM_(0),
    HLT_Photon75_R9Id90_HE10_IsoM_(0),
    HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_(0),
    HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_(0),
    HLT_Photon90_R9Id90_HE10_IsoM_(0),
    HLT_Photon120_R9Id90_HE10_IsoM_(0),
    HLT_Photon165_R9Id90_HE10_IsoM_(0),
    HLT_Photon90_CaloIdL_PFHT700_(0),
    HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_(0),
    HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_(0),
    HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_(0),
    HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_(0),
    HLT_Photon35_TwoProngs35_(0),
    HLT_IsoMu24_TwoProngs35_(0),
    HLT_Dimuon0_Jpsi_L1_NoOS_(0),
    HLT_Dimuon0_Jpsi_NoVertexing_NoOS_(0),
    HLT_Dimuon0_Jpsi_(0),
    HLT_Dimuon0_Jpsi_NoVertexing_(0),
    HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_(0),
    HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_(0),
    HLT_Dimuon0_Jpsi3p5_Muon2_(0),
    HLT_Dimuon0_Upsilon_L1_4p5_(0),
    HLT_Dimuon0_Upsilon_L1_5_(0),
    HLT_Dimuon0_Upsilon_L1_4p5NoOS_(0),
    HLT_Dimuon0_Upsilon_L1_4p5er2p0_(0),
    HLT_Dimuon0_Upsilon_L1_4p5er2p0M_(0),
    HLT_Dimuon0_Upsilon_NoVertexing_(0),
    HLT_Dimuon0_Upsilon_L1_5M_(0),
    HLT_Dimuon0_LowMass_L1_0er1p5R_(0),
    HLT_Dimuon0_LowMass_L1_0er1p5_(0),
    HLT_Dimuon0_LowMass_(0),
    HLT_Dimuon0_LowMass_L1_4_(0),
    HLT_Dimuon0_LowMass_L1_4R_(0),
    HLT_Dimuon0_LowMass_L1_TM530_(0),
    HLT_Dimuon0_Upsilon_Muon_L1_TM0_(0),
    HLT_Dimuon0_Upsilon_Muon_NoL1Mass_(0),
    HLT_TripleMu_5_3_3_Mass3p8_DZ_(0),
    HLT_TripleMu_10_5_5_DZ_(0),
    HLT_TripleMu_12_10_5_(0),
    HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_(0),
    HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_(0),
    HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_(0),
    HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_(0),
    HLT_DoubleMu3_DZ_PFMET50_PFMHT60_(0),
    HLT_DoubleMu3_DZ_PFMET70_PFMHT70_(0),
    HLT_DoubleMu3_DZ_PFMET90_PFMHT90_(0),
    HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_(0),
    HLT_DoubleMu4_Jpsi_Displaced_(0),
    HLT_DoubleMu4_Jpsi_NoVertexing_(0),
    HLT_DoubleMu4_JpsiTrkTrk_Displaced_(0),
    HLT_DoubleMu43NoFiltersNoVtx_(0),
    HLT_DoubleMu48NoFiltersNoVtx_(0),
    HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_(0),
    HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_(0),
    HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_(0),
    HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_(0),
    HLT_DoubleMu33NoFiltersNoVtxDisplaced_(0),
    HLT_DoubleMu40NoFiltersNoVtxDisplaced_(0),
    HLT_DoubleMu20_7_Mass0to30_L1_DM4_(0),
    HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_(0),
    HLT_HT425_(0),
    HLT_HT430_DisplacedDijet40_DisplacedTrack_(0),
    HLT_HT500_DisplacedDijet40_DisplacedTrack_(0),
    HLT_HT430_DisplacedDijet60_DisplacedTrack_(0),
    HLT_HT400_DisplacedDijet40_DisplacedTrack_(0),
    HLT_HT650_DisplacedDijet60_Inclusive_(0),
    HLT_HT550_DisplacedDijet60_Inclusive_(0),
    HLT_DiJet110_35_Mjj650_PFMET110_(0),
    HLT_DiJet110_35_Mjj650_PFMET120_(0),
    HLT_DiJet110_35_Mjj650_PFMET130_(0),
    HLT_TripleJet110_35_35_Mjj650_PFMET110_(0),
    HLT_TripleJet110_35_35_Mjj650_PFMET120_(0),
    HLT_TripleJet110_35_35_Mjj650_PFMET130_(0),
    HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_(0),
    HLT_Ele28_eta2p1_WPTight_Gsf_HT150_(0),
    HLT_Ele28_HighEta_SC20_Mass55_(0),
    HLT_DoubleMu20_7_Mass0to30_Photon23_(0),
    HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_(0),
    HLT_Ele15_IsoVVVL_PFHT450_PFMET50_(0),
    HLT_Ele15_IsoVVVL_PFHT450_(0),
    HLT_Ele50_IsoVVVL_PFHT450_(0),
    HLT_Ele15_IsoVVVL_PFHT600_(0),
    HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_(0),
    HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_(0),
    HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_(0),
    HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_(0),
    HLT_Mu15_IsoVVVL_PFHT450_PFMET50_(0),
    HLT_Mu15_IsoVVVL_PFHT450_(0),
    HLT_Mu50_IsoVVVL_PFHT450_(0),
    HLT_Mu15_IsoVVVL_PFHT600_(0),
    HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_(0),
    HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_(0),
    HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_(0),
    HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_(0),
    HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_(0),
    HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_(0),
    HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_(0),
    HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_(0),
    HLT_Dimuon10_PsiPrime_Barrel_Seagulls_(0),
    HLT_Dimuon20_Jpsi_Barrel_Seagulls_(0),
    HLT_Dimuon12_Upsilon_y1p4_(0),
    HLT_Dimuon14_Phi_Barrel_Seagulls_(0),
    HLT_Dimuon18_PsiPrime_(0),
    HLT_Dimuon25_Jpsi_(0),
    HLT_Dimuon18_PsiPrime_noCorrL1_(0),
    HLT_Dimuon24_Upsilon_noCorrL1_(0),
    HLT_Dimuon24_Phi_noCorrL1_(0),
    HLT_Dimuon25_Jpsi_noCorrL1_(0),
    HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_(0),
    HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_(0),
    HLT_DiMu9_Ele9_CaloIdL_TrackIdL_(0),
    HLT_DoubleIsoMu20_eta2p1_(0),
    HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_(0),
    HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_(0),
    HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_(0),
    HLT_Mu8_(0),
    HLT_Mu17_(0),
    HLT_Mu19_(0),
    HLT_Mu17_Photon30_IsoCaloId_(0),
    HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_(0),
    HLT_Ele8_CaloIdM_TrackIdM_PFJet30_(0),
    HLT_Ele17_CaloIdM_TrackIdM_PFJet30_(0),
    HLT_Ele23_CaloIdM_TrackIdM_PFJet30_(0),
    HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_(0),
    HLT_Ele115_CaloIdVT_GsfTrkIdT_(0),
    HLT_Ele135_CaloIdVT_GsfTrkIdT_(0),
    HLT_Ele145_CaloIdVT_GsfTrkIdT_(0),
    HLT_Ele200_CaloIdVT_GsfTrkIdT_(0),
    HLT_Ele250_CaloIdVT_GsfTrkIdT_(0),
    HLT_Ele300_CaloIdVT_GsfTrkIdT_(0),
    HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_(0),
    HLT_PFHT330PT30_QuadPFJet_75_60_45_40_(0),
    HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_(0),
    HLT_PFHT400_SixPFJet32_(0),
    HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_(0),
    HLT_PFHT450_SixPFJet36_(0),
    HLT_PFHT350_(0),
    HLT_PFHT350MinPFJet15_(0),
    HLT_Photon60_R9Id90_CaloIdL_IsoL_(0),
    HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_(0),
    HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_(0),
    HLT_ECALHT800_(0),
    HLT_DiSC30_18_EIso_AND_HE_Mass70_(0),
    HLT_Physics_(0),
    HLT_Physics_part0_(0),
    HLT_Physics_part1_(0),
    HLT_Physics_part2_(0),
    HLT_Physics_part3_(0),
    HLT_Physics_part4_(0),
    HLT_Physics_part5_(0),
    HLT_Physics_part6_(0),
    HLT_Physics_part7_(0),
    HLT_Random_(0),
    HLT_ZeroBias_(0),
    HLT_ZeroBias_Alignment_(0),
    HLT_ZeroBias_part0_(0),
    HLT_ZeroBias_part1_(0),
    HLT_ZeroBias_part2_(0),
    HLT_ZeroBias_part3_(0),
    HLT_ZeroBias_part4_(0),
    HLT_ZeroBias_part5_(0),
    HLT_ZeroBias_part6_(0),
    HLT_ZeroBias_part7_(0),
    HLT_AK4CaloJet30_(0),
    HLT_AK4CaloJet40_(0),
    HLT_AK4CaloJet50_(0),
    HLT_AK4CaloJet80_(0),
    HLT_AK4CaloJet100_(0),
    HLT_AK4CaloJet120_(0),
    HLT_AK4PFJet30_(0),
    HLT_AK4PFJet50_(0),
    HLT_AK4PFJet80_(0),
    HLT_AK4PFJet100_(0),
    HLT_AK4PFJet120_(0),
    HLT_SinglePhoton10_Eta3p1ForPPRef_(0),
    HLT_SinglePhoton20_Eta3p1ForPPRef_(0),
    HLT_SinglePhoton30_Eta3p1ForPPRef_(0),
    HLT_Photon20_HoverELoose_(0),
    HLT_Photon30_HoverELoose_(0),
    HLT_EcalCalibration_(0),
    HLT_HcalCalibration_(0),
    HLT_L1UnpairedBunchBptxMinus_(0),
    HLT_L1UnpairedBunchBptxPlus_(0),
    HLT_L1NotBptxOR_(0),
    HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_(0),
    HLT_CDC_L2cosmic_5_er1p0_(0),
    HLT_CDC_L2cosmic_5p5_er1p0_(0),
    HLT_HcalNZS_(0),
    HLT_HcalPhiSym_(0),
    HLT_HcalIsolatedbunch_(0),
    HLT_IsoTrackHB_(0),
    HLT_IsoTrackHE_(0),
    HLT_ZeroBias_FirstCollisionAfterAbortGap_(0),
    HLT_ZeroBias_IsolatedBunches_(0),
    HLT_ZeroBias_FirstCollisionInTrain_(0),
    HLT_ZeroBias_LastCollisionInTrain_(0),
    HLT_ZeroBias_FirstBXAfterTrain_(0),
    HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_(0),
    HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_(0),
    HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_(0),
    HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_(0),
    HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_(0),
    HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_(0),
    HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_(0),
    HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_(0),
    HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_(0),
    HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_(0),
    HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_(0),
    HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_(0),
    HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_(0),
    HLT_Rsq0p35_(0),
    HLT_Rsq0p40_(0),
    HLT_RsqMR300_Rsq0p09_MR200_(0),
    HLT_RsqMR320_Rsq0p09_MR200_(0),
    HLT_RsqMR300_Rsq0p09_MR200_4jet_(0),
    HLT_RsqMR320_Rsq0p09_MR200_4jet_(0),
    HLT_IsoMu27_MET90_(0),
    HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_(0),
    HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_(0),
    HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_(0),
    HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_(0),
    HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_(0),
    HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_(0),
    HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_(0),
    HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_(0),
    HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_(0),
    HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_(0),
    HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_(0),
    HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_(0),
    HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_(0),
    HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_(0),
    HLT_PFMET100_PFMHT100_IDTight_PFHT60_(0),
    HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_(0),
    HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_(0),
    HLT_Mu18_Mu9_SameSign_(0),
    HLT_Mu18_Mu9_SameSign_DZ_(0),
    HLT_Mu18_Mu9_(0),
    HLT_Mu18_Mu9_DZ_(0),
    HLT_Mu20_Mu10_SameSign_(0),
    HLT_Mu20_Mu10_SameSign_DZ_(0),
    HLT_Mu20_Mu10_(0),
    HLT_Mu20_Mu10_DZ_(0),
    HLT_Mu23_Mu12_SameSign_(0),
    HLT_Mu23_Mu12_SameSign_DZ_(0),
    HLT_Mu23_Mu12_(0),
    HLT_Mu23_Mu12_DZ_(0),
    HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_(0),
    HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_(0),
    HLT_DoubleMu3_DCA_PFMET50_PFMHT60_(0),
    HLT_TripleMu_5_3_3_Mass3p8_DCA_(0),
    HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_(0),
    HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_(0),
    HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_(0),
    HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_(0),
    HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_(0),
    HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_(0),
    HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_(0),
    HLT_QuadPFJet98_83_71_15_(0),
    HLT_QuadPFJet103_88_75_15_(0),
    HLT_QuadPFJet105_88_76_15_(0),
    HLT_QuadPFJet111_90_80_15_(0),
    HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_(0),
    HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_(0),
    HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_(0),
    HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_(0),
    HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_(0),
    HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_(0),
    HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_(0),
    HLT_Mu12_IP6_part0_(0),
    HLT_Mu12_IP6_part1_(0),
    HLT_Mu12_IP6_part2_(0),
    HLT_Mu12_IP6_part3_(0),
    HLT_Mu12_IP6_part4_(0),
    HLT_Mu9_IP5_part0_(0),
    HLT_Mu9_IP5_part1_(0),
    HLT_Mu9_IP5_part2_(0),
    HLT_Mu9_IP5_part3_(0),
    HLT_Mu9_IP5_part4_(0),
    HLT_Mu7_IP4_part0_(0),
    HLT_Mu7_IP4_part1_(0),
    HLT_Mu7_IP4_part2_(0),
    HLT_Mu7_IP4_part3_(0),
    HLT_Mu7_IP4_part4_(0),
    HLT_Mu9_IP4_part0_(0),
    HLT_Mu9_IP4_part1_(0),
    HLT_Mu9_IP4_part2_(0),
    HLT_Mu9_IP4_part3_(0),
    HLT_Mu9_IP4_part4_(0),
    HLT_Mu8_IP5_part0_(0),
    HLT_Mu8_IP5_part1_(0),
    HLT_Mu8_IP5_part2_(0),
    HLT_Mu8_IP5_part3_(0),
    HLT_Mu8_IP5_part4_(0),
    HLT_Mu8_IP6_part0_(0),
    HLT_Mu8_IP6_part1_(0),
    HLT_Mu8_IP6_part2_(0),
    HLT_Mu8_IP6_part3_(0),
    HLT_Mu8_IP6_part4_(0),
    HLT_Mu9_IP6_part0_(0),
    HLT_Mu9_IP6_part1_(0),
    HLT_Mu9_IP6_part2_(0),
    HLT_Mu9_IP6_part3_(0),
    HLT_Mu9_IP6_part4_(0),
    HLT_Mu8_IP3_part0_(0),
    HLT_Mu8_IP3_part1_(0),
    HLT_Mu8_IP3_part2_(0),
    HLT_Mu8_IP3_part3_(0),
    HLT_Mu8_IP3_part4_(0),
    HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_(0),
    HLT_TrkMu6NoFiltersNoVtx_(0),
    HLT_TrkMu16NoFiltersNoVtx_(0),
    HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_(0),
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
    Flag_ecalBadCalibFilter_(0),
    Flag_goodVertices_(0),
    Flag_eeBadScFilter_(0),
    Flag_ecalLaserCorrFilter_(0),
    Flag_trkPOGFilters_(0),
    Flag_chargedHadronTrackResolutionFilter_(0),
    Flag_muonBadTrackFilter_(0),
    Flag_BadChargedCandidateFilter_(0),
    Flag_BadPFMuonFilter_(0),
    Flag_BadChargedCandidateSummer16Filter_(0),
    Flag_BadPFMuonSummer16Filter_(0),
    Flag_trkPOG_manystripclus53X_(0),
    Flag_trkPOG_toomanystripclus53X_(0),
    Flag_trkPOG_logErrorTooManyClusters_(0),
    Flag_METFilters_(0),
    L1Reco_step_(0),
    are_Jet_loaded_(0), Jet_(),
    are_IsoTrack_loaded_(0), IsoTrack_(),
    are_GenJetAK8_loaded_(0), GenJetAK8_(),
    are_GenVisTau_loaded_(0), GenVisTau_(),
    are_CaloMET_loaded_(0), CaloMET_(),
    are_GenDressedLepton_loaded_(0), GenDressedLepton_(),
    are_LHEPart_loaded_(0), LHEPart_(),
    are_PV_loaded_(0), PV_(),
    are_Generator_loaded_(0), Generator_(),
    are_TrigObj_loaded_(0), TrigObj_(),
    are_Photon_loaded_(0), Photon_(),
    are_GenJet_loaded_(0), GenJet_(),
    are_RawMET_loaded_(0), RawMET_(),
    are_btagWeight_loaded_(0), btagWeight_(),
    are_SoftActivityJet_loaded_(0), SoftActivityJet_(),
    are_L1simulation_loaded_(0), L1simulation_(),
    are_GenPart_loaded_(0), GenPart_(),
    are_LHE_loaded_(0), LHE_(),
    are_TkMET_loaded_(0), TkMET_(),
    are_Tau_loaded_(0), Tau_(),
    are_Flag_loaded_(0), Flag_(),
    are_PuppiMET_loaded_(0), PuppiMET_(),
    are_Muon_loaded_(0), Muon_(),
    are_L1Reco_loaded_(0), L1Reco_(),
    are_OtherPV_loaded_(0), OtherPV_(),
    are_HLT_loaded_(0), HLT_(),
    are_SV_loaded_(0), SV_(),
    are_Electron_loaded_(0), Electron_(),
    are_MET_loaded_(0), MET_(),
    are_GenMET_loaded_(0), GenMET_(),
    are_FatJet_loaded_(0), FatJet_(),
    are_SubJet_loaded_(0), SubJet_(),
    are_SubGenJetAK8_loaded_(0), SubGenJetAK8_(),
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
    tree_->SetBranchStatus("nSubGenJetAK8", 1); tree_->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8);
    tree_->SetBranchStatus("nGenVisTau", 1); tree_->SetBranchAddress("nGenVisTau", &nGenVisTau);
    tree_->SetBranchStatus("genWeight", 1); tree_->SetBranchAddress("genWeight", &genWeight);
    tree_->SetBranchStatus("nLHEPdfWeight", 1); tree_->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight);
    tree_->SetBranchStatus("LHEPdfWeight", 1); tree_->SetBranchAddress("LHEPdfWeight", &LHEPdfWeight);
    tree_->SetBranchStatus("nLHEScaleWeight", 1); tree_->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight);
    tree_->SetBranchStatus("LHEScaleWeight", 1); tree_->SetBranchAddress("LHEScaleWeight", &LHEScaleWeight);
    tree_->SetBranchStatus("nPSWeight", 1); tree_->SetBranchAddress("nPSWeight", &nPSWeight);
    tree_->SetBranchStatus("PSWeight", 1); tree_->SetBranchAddress("PSWeight", &PSWeight);
    tree_->SetBranchStatus("nIsoTrack", 1); tree_->SetBranchAddress("nIsoTrack", &nIsoTrack);
    tree_->SetBranchStatus("nJet", 1); tree_->SetBranchAddress("nJet", &nJet);
    tree_->SetBranchStatus("nLHEPart", 1); tree_->SetBranchAddress("nLHEPart", &nLHEPart);
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
    IsoTrack_.clear();
    GenJetAK8_.clear();
    GenVisTau_.clear();
    
    GenDressedLepton_.clear();
    LHEPart_.clear();
    
    
    TrigObj_.clear();
    Photon_.clear();
    GenJet_.clear();
    
    
    
    
    GenPart_.clear();
    
    
    Tau_.clear();
    
    
    Muon_.clear();
    
    OtherPV_.clear();
    
    SV_.clear();
    Electron_.clear();
    
    
    FatJet_.clear();
    SubJet_.clear();
    SubGenJetAK8_.clear();
    
    
    current_entry_++;
    if(current_entry_ < entries_){
      tree_->GetEntry(current_entry_);
      return true;
    }
    return false;
  }

  void loadJets(){
    if(!are_Jet_loaded_){
      tree_->SetBranchStatus("Jet_area", 1); tree_->SetBranchAddress("Jet_area", &Jet_area_);
      tree_->SetBranchStatus("Jet_btagCMVA", 1); tree_->SetBranchAddress("Jet_btagCMVA", &Jet_btagCMVA_);
      tree_->SetBranchStatus("Jet_btagCSVV2", 1); tree_->SetBranchAddress("Jet_btagCSVV2", &Jet_btagCSVV2_);
      tree_->SetBranchStatus("Jet_btagDeepB", 1); tree_->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB_);
      tree_->SetBranchStatus("Jet_btagDeepC", 1); tree_->SetBranchAddress("Jet_btagDeepC", &Jet_btagDeepC_);
      tree_->SetBranchStatus("Jet_btagDeepFlavB", 1); tree_->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB_);
      tree_->SetBranchStatus("Jet_chEmEF", 1); tree_->SetBranchAddress("Jet_chEmEF", &Jet_chEmEF_);
      tree_->SetBranchStatus("Jet_chHEF", 1); tree_->SetBranchAddress("Jet_chHEF", &Jet_chHEF_);
      tree_->SetBranchStatus("Jet_eta", 1); tree_->SetBranchAddress("Jet_eta", &Jet_eta_);
      tree_->SetBranchStatus("Jet_mass", 1); tree_->SetBranchAddress("Jet_mass", &Jet_mass_);
      tree_->SetBranchStatus("Jet_muEF", 1); tree_->SetBranchAddress("Jet_muEF", &Jet_muEF_);
      tree_->SetBranchStatus("Jet_neEmEF", 1); tree_->SetBranchAddress("Jet_neEmEF", &Jet_neEmEF_);
      tree_->SetBranchStatus("Jet_neHEF", 1); tree_->SetBranchAddress("Jet_neHEF", &Jet_neHEF_);
      tree_->SetBranchStatus("Jet_phi", 1); tree_->SetBranchAddress("Jet_phi", &Jet_phi_);
      tree_->SetBranchStatus("Jet_pt", 1); tree_->SetBranchAddress("Jet_pt", &Jet_pt_);
      tree_->SetBranchStatus("Jet_qgl", 1); tree_->SetBranchAddress("Jet_qgl", &Jet_qgl_);
      tree_->SetBranchStatus("Jet_rawFactor", 1); tree_->SetBranchAddress("Jet_rawFactor", &Jet_rawFactor_);
      tree_->SetBranchStatus("Jet_bRegCorr", 1); tree_->SetBranchAddress("Jet_bRegCorr", &Jet_bRegCorr_);
      tree_->SetBranchStatus("Jet_bRegRes", 1); tree_->SetBranchAddress("Jet_bRegRes", &Jet_bRegRes_);
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
  
  void loadIsotracks(){
    if(!are_IsoTrack_loaded_){
      tree_->SetBranchStatus("IsoTrack_dxy", 1); tree_->SetBranchAddress("IsoTrack_dxy", &IsoTrack_dxy_);
      tree_->SetBranchStatus("IsoTrack_dz", 1); tree_->SetBranchAddress("IsoTrack_dz", &IsoTrack_dz_);
      tree_->SetBranchStatus("IsoTrack_eta", 1); tree_->SetBranchAddress("IsoTrack_eta", &IsoTrack_eta_);
      tree_->SetBranchStatus("IsoTrack_pfRelIso03_all", 1); tree_->SetBranchAddress("IsoTrack_pfRelIso03_all", &IsoTrack_pfRelIso03_all_);
      tree_->SetBranchStatus("IsoTrack_pfRelIso03_chg", 1); tree_->SetBranchAddress("IsoTrack_pfRelIso03_chg", &IsoTrack_pfRelIso03_chg_);
      tree_->SetBranchStatus("IsoTrack_phi", 1); tree_->SetBranchAddress("IsoTrack_phi", &IsoTrack_phi_);
      tree_->SetBranchStatus("IsoTrack_pt", 1); tree_->SetBranchAddress("IsoTrack_pt", &IsoTrack_pt_);
      tree_->SetBranchStatus("IsoTrack_miniPFRelIso_all", 1); tree_->SetBranchAddress("IsoTrack_miniPFRelIso_all", &IsoTrack_miniPFRelIso_all_);
      tree_->SetBranchStatus("IsoTrack_miniPFRelIso_chg", 1); tree_->SetBranchAddress("IsoTrack_miniPFRelIso_chg", &IsoTrack_miniPFRelIso_chg_);
      tree_->SetBranchStatus("IsoTrack_pdgId", 1); tree_->SetBranchAddress("IsoTrack_pdgId", &IsoTrack_pdgId_);
      tree_->SetBranchStatus("IsoTrack_isHighPurityTrack", 1); tree_->SetBranchAddress("IsoTrack_isHighPurityTrack", &IsoTrack_isHighPurityTrack_);
      tree_->SetBranchStatus("IsoTrack_isPFcand", 1); tree_->SetBranchAddress("IsoTrack_isPFcand", &IsoTrack_isPFcand_);
      are_IsoTrack_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadGenjetak8s(){
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
  
  void loadGenvistaus(){
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
  
  void loadGendressedleptons(){
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
  
  void loadLheparts(){
    if(!are_LHEPart_loaded_){
      tree_->SetBranchStatus("LHEPart_pt", 1); tree_->SetBranchAddress("LHEPart_pt", &LHEPart_pt_);
      tree_->SetBranchStatus("LHEPart_eta", 1); tree_->SetBranchAddress("LHEPart_eta", &LHEPart_eta_);
      tree_->SetBranchStatus("LHEPart_phi", 1); tree_->SetBranchAddress("LHEPart_phi", &LHEPart_phi_);
      tree_->SetBranchStatus("LHEPart_mass", 1); tree_->SetBranchAddress("LHEPart_mass", &LHEPart_mass_);
      tree_->SetBranchStatus("LHEPart_pdgId", 1); tree_->SetBranchAddress("LHEPart_pdgId", &LHEPart_pdgId_);
      are_LHEPart_loaded_ = true;
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
  
  void loadTrigobjs(){
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
  
  void loadPhotons(){
    if(!are_Photon_loaded_){
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
      tree_->SetBranchStatus("Photon_cutBasedBitmap", 1); tree_->SetBranchAddress("Photon_cutBasedBitmap", &Photon_cutBasedBitmap_);
      tree_->SetBranchStatus("Photon_electronIdx", 1); tree_->SetBranchAddress("Photon_electronIdx", &Photon_electronIdx_);
      tree_->SetBranchStatus("Photon_jetIdx", 1); tree_->SetBranchAddress("Photon_jetIdx", &Photon_jetIdx_);
      tree_->SetBranchStatus("Photon_pdgId", 1); tree_->SetBranchAddress("Photon_pdgId", &Photon_pdgId_);
      tree_->SetBranchStatus("Photon_vidNestedWPBitmap", 1); tree_->SetBranchAddress("Photon_vidNestedWPBitmap", &Photon_vidNestedWPBitmap_);
      tree_->SetBranchStatus("Photon_electronVeto", 1); tree_->SetBranchAddress("Photon_electronVeto", &Photon_electronVeto_);
      tree_->SetBranchStatus("Photon_isScEtaEB", 1); tree_->SetBranchAddress("Photon_isScEtaEB", &Photon_isScEtaEB_);
      tree_->SetBranchStatus("Photon_isScEtaEE", 1); tree_->SetBranchAddress("Photon_isScEtaEE", &Photon_isScEtaEE_);
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
  
  void loadGenjets(){
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
  
  void loadBtagweight(){
    if(!are_btagWeight_loaded_){
      tree_->SetBranchStatus("btagWeight_CSVV2", 1); tree_->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2_);
      tree_->SetBranchStatus("btagWeight_DeepCSVB", 1); tree_->SetBranchAddress("btagWeight_DeepCSVB", &btagWeight_DeepCSVB_);
      are_btagWeight_loaded_ = true;
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
  
  void loadGenparts(){
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
      tree_->SetBranchStatus("LHEPart_pt", 1); tree_->SetBranchAddress("LHEPart_pt", &LHEPart_pt_);
      tree_->SetBranchStatus("LHEPart_eta", 1); tree_->SetBranchAddress("LHEPart_eta", &LHEPart_eta_);
      tree_->SetBranchStatus("LHEPart_phi", 1); tree_->SetBranchAddress("LHEPart_phi", &LHEPart_phi_);
      tree_->SetBranchStatus("LHEPart_mass", 1); tree_->SetBranchAddress("LHEPart_mass", &LHEPart_mass_);
      tree_->SetBranchStatus("LHEPart_pdgId", 1); tree_->SetBranchAddress("LHEPart_pdgId", &LHEPart_pdgId_);
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
  
  void loadTaus(){
    if(!are_Tau_loaded_){
      tree_->SetBranchStatus("Tau_chargedIso", 1); tree_->SetBranchAddress("Tau_chargedIso", &Tau_chargedIso_);
      tree_->SetBranchStatus("Tau_dxy", 1); tree_->SetBranchAddress("Tau_dxy", &Tau_dxy_);
      tree_->SetBranchStatus("Tau_dz", 1); tree_->SetBranchAddress("Tau_dz", &Tau_dz_);
      tree_->SetBranchStatus("Tau_eta", 1); tree_->SetBranchAddress("Tau_eta", &Tau_eta_);
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
      tree_->SetBranchStatus("Tau_rawIsodR03", 1); tree_->SetBranchAddress("Tau_rawIsodR03", &Tau_rawIsodR03_);
      tree_->SetBranchStatus("Tau_charge", 1); tree_->SetBranchAddress("Tau_charge", &Tau_charge_);
      tree_->SetBranchStatus("Tau_decayMode", 1); tree_->SetBranchAddress("Tau_decayMode", &Tau_decayMode_);
      tree_->SetBranchStatus("Tau_jetIdx", 1); tree_->SetBranchAddress("Tau_jetIdx", &Tau_jetIdx_);
      tree_->SetBranchStatus("Tau_rawAntiEleCat", 1); tree_->SetBranchAddress("Tau_rawAntiEleCat", &Tau_rawAntiEleCat_);
      tree_->SetBranchStatus("Tau_idAntiEle", 1); tree_->SetBranchAddress("Tau_idAntiEle", &Tau_idAntiEle_);
      tree_->SetBranchStatus("Tau_idAntiMu", 1); tree_->SetBranchAddress("Tau_idAntiMu", &Tau_idAntiMu_);
      tree_->SetBranchStatus("Tau_idDecayMode", 1); tree_->SetBranchAddress("Tau_idDecayMode", &Tau_idDecayMode_);
      tree_->SetBranchStatus("Tau_idDecayModeNewDMs", 1); tree_->SetBranchAddress("Tau_idDecayModeNewDMs", &Tau_idDecayModeNewDMs_);
      tree_->SetBranchStatus("Tau_cleanmask", 1); tree_->SetBranchAddress("Tau_cleanmask", &Tau_cleanmask_);
      tree_->SetBranchStatus("Tau_genPartIdx", 1); tree_->SetBranchAddress("Tau_genPartIdx", &Tau_genPartIdx_);
      tree_->SetBranchStatus("Tau_genPartFlav", 1); tree_->SetBranchAddress("Tau_genPartFlav", &Tau_genPartFlav_);
      are_Tau_loaded_ = true;
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
      tree_->SetBranchStatus("Flag_ecalBadCalibFilter", 1); tree_->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter_);
      tree_->SetBranchStatus("Flag_goodVertices", 1); tree_->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices_);
      tree_->SetBranchStatus("Flag_eeBadScFilter", 1); tree_->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter_);
      tree_->SetBranchStatus("Flag_ecalLaserCorrFilter", 1); tree_->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter_);
      tree_->SetBranchStatus("Flag_trkPOGFilters", 1); tree_->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters_);
      tree_->SetBranchStatus("Flag_chargedHadronTrackResolutionFilter", 1); tree_->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter_);
      tree_->SetBranchStatus("Flag_muonBadTrackFilter", 1); tree_->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter_);
      tree_->SetBranchStatus("Flag_BadChargedCandidateFilter", 1); tree_->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter_);
      tree_->SetBranchStatus("Flag_BadPFMuonFilter", 1); tree_->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter_);
      tree_->SetBranchStatus("Flag_BadChargedCandidateSummer16Filter", 1); tree_->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter_);
      tree_->SetBranchStatus("Flag_BadPFMuonSummer16Filter", 1); tree_->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter_);
      tree_->SetBranchStatus("Flag_trkPOG_manystripclus53X", 1); tree_->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X_);
      tree_->SetBranchStatus("Flag_trkPOG_toomanystripclus53X", 1); tree_->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X_);
      tree_->SetBranchStatus("Flag_trkPOG_logErrorTooManyClusters", 1); tree_->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters_);
      tree_->SetBranchStatus("Flag_METFilters", 1); tree_->SetBranchAddress("Flag_METFilters", &Flag_METFilters_);
      are_Flag_loaded_ = true;
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
  
  void loadMuons(){
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
  
  void loadL1Reco(){
    if(!are_L1Reco_loaded_){
      tree_->SetBranchStatus("L1Reco_step", 1); tree_->SetBranchAddress("L1Reco_step", &L1Reco_step_);
      are_L1Reco_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadOtherpvs(){
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
      tree_->SetBranchStatus("HLT_AK8PFJet380_TrimMass30", 1); tree_->SetBranchAddress("HLT_AK8PFJet380_TrimMass30", &HLT_AK8PFJet380_TrimMass30_);
      tree_->SetBranchStatus("HLT_AK8PFJet400_TrimMass30", 1); tree_->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30_);
      tree_->SetBranchStatus("HLT_AK8PFJet420_TrimMass30", 1); tree_->SetBranchAddress("HLT_AK8PFJet420_TrimMass30", &HLT_AK8PFJet420_TrimMass30_);
      tree_->SetBranchStatus("HLT_AK8PFHT750_TrimMass50", 1); tree_->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50_);
      tree_->SetBranchStatus("HLT_AK8PFHT800_TrimMass50", 1); tree_->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50_);
      tree_->SetBranchStatus("HLT_AK8PFHT850_TrimMass50", 1); tree_->SetBranchAddress("HLT_AK8PFHT850_TrimMass50", &HLT_AK8PFHT850_TrimMass50_);
      tree_->SetBranchStatus("HLT_AK8PFHT900_TrimMass50", 1); tree_->SetBranchAddress("HLT_AK8PFHT900_TrimMass50", &HLT_AK8PFHT900_TrimMass50_);
      tree_->SetBranchStatus("HLT_CaloJet500_NoJetID", 1); tree_->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID_);
      tree_->SetBranchStatus("HLT_CaloJet550_NoJetID", 1); tree_->SetBranchAddress("HLT_CaloJet550_NoJetID", &HLT_CaloJet550_NoJetID_);
      tree_->SetBranchStatus("HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL", 1); tree_->SetBranchAddress("HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL", &HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_);
      tree_->SetBranchStatus("HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon", 1); tree_->SetBranchAddress("HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon", &HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_);
      tree_->SetBranchStatus("HLT_Trimuon5_3p5_2_Upsilon_Muon", 1); tree_->SetBranchAddress("HLT_Trimuon5_3p5_2_Upsilon_Muon", &HLT_Trimuon5_3p5_2_Upsilon_Muon_);
      tree_->SetBranchStatus("HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon", 1); tree_->SetBranchAddress("HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon", &HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_);
      tree_->SetBranchStatus("HLT_DoubleEle25_CaloIdL_MW", 1); tree_->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MW", &HLT_DoubleEle25_CaloIdL_MW_);
      tree_->SetBranchStatus("HLT_DoubleEle27_CaloIdL_MW", 1); tree_->SetBranchAddress("HLT_DoubleEle27_CaloIdL_MW", &HLT_DoubleEle27_CaloIdL_MW_);
      tree_->SetBranchStatus("HLT_DoubleEle33_CaloIdL_MW", 1); tree_->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW_);
      tree_->SetBranchStatus("HLT_DoubleEle24_eta2p1_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_DoubleEle24_eta2p1_WPTight_Gsf", &HLT_DoubleEle24_eta2p1_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", 1); tree_->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_);
      tree_->SetBranchStatus("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", 1); tree_->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_);
      tree_->SetBranchStatus("HLT_Ele27_Ele37_CaloIdL_MW", 1); tree_->SetBranchAddress("HLT_Ele27_Ele37_CaloIdL_MW", &HLT_Ele27_Ele37_CaloIdL_MW_);
      tree_->SetBranchStatus("HLT_Mu27_Ele37_CaloIdL_MW", 1); tree_->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_MW", &HLT_Mu27_Ele37_CaloIdL_MW_);
      tree_->SetBranchStatus("HLT_Mu37_Ele27_CaloIdL_MW", 1); tree_->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_MW", &HLT_Mu37_Ele27_CaloIdL_MW_);
      tree_->SetBranchStatus("HLT_Mu37_TkMu27", 1); tree_->SetBranchAddress("HLT_Mu37_TkMu27", &HLT_Mu37_TkMu27_);
      tree_->SetBranchStatus("HLT_DoubleMu4_3_Bs", 1); tree_->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs_);
      tree_->SetBranchStatus("HLT_DoubleMu4_3_Jpsi", 1); tree_->SetBranchAddress("HLT_DoubleMu4_3_Jpsi", &HLT_DoubleMu4_3_Jpsi_);
      tree_->SetBranchStatus("HLT_DoubleMu4_JpsiTrk_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced_);
      tree_->SetBranchStatus("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_);
      tree_->SetBranchStatus("HLT_DoubleMu3_Trk_Tau3mu", 1); tree_->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu_);
      tree_->SetBranchStatus("HLT_DoubleMu3_TkMu_DsTau3Mu", 1); tree_->SetBranchAddress("HLT_DoubleMu3_TkMu_DsTau3Mu", &HLT_DoubleMu3_TkMu_DsTau3Mu_);
      tree_->SetBranchStatus("HLT_DoubleMu4_PsiPrimeTrk_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced_);
      tree_->SetBranchStatus("HLT_DoubleMu4_Mass3p8_DZ_PFHT350", 1); tree_->SetBranchAddress("HLT_DoubleMu4_Mass3p8_DZ_PFHT350", &HLT_DoubleMu4_Mass3p8_DZ_PFHT350_);
      tree_->SetBranchStatus("HLT_Mu3_PFJet40", 1); tree_->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40_);
      tree_->SetBranchStatus("HLT_Mu7p5_L2Mu2_Jpsi", 1); tree_->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi_);
      tree_->SetBranchStatus("HLT_Mu7p5_L2Mu2_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track2_Jpsi", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track3p5_Jpsi", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track7_Jpsi", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track2_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track3p5_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon_);
      tree_->SetBranchStatus("HLT_Mu7p5_Track7_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon_);
      tree_->SetBranchStatus("HLT_Mu3_L1SingleMu5orSingleMu7", 1); tree_->SetBranchAddress("HLT_Mu3_L1SingleMu5orSingleMu7", &HLT_Mu3_L1SingleMu5orSingleMu7_);
      tree_->SetBranchStatus("HLT_DoublePhoton33_CaloIdL", 1); tree_->SetBranchAddress("HLT_DoublePhoton33_CaloIdL", &HLT_DoublePhoton33_CaloIdL_);
      tree_->SetBranchStatus("HLT_DoublePhoton70", 1); tree_->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70_);
      tree_->SetBranchStatus("HLT_DoublePhoton85", 1); tree_->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85_);
      tree_->SetBranchStatus("HLT_Ele20_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele20_WPTight_Gsf", &HLT_Ele20_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele15_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele15_WPLoose_Gsf", &HLT_Ele15_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele17_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele17_WPLoose_Gsf", &HLT_Ele17_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele20_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele20_WPLoose_Gsf", &HLT_Ele20_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_Ele20_eta2p1_WPLoose_Gsf", 1); tree_->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf", &HLT_Ele20_eta2p1_WPLoose_Gsf_);
      tree_->SetBranchStatus("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", 1); tree_->SetBranchAddress("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_);
      tree_->SetBranchStatus("HLT_Ele27_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele28_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele28_WPTight_Gsf", &HLT_Ele28_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele30_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele35_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele35_WPTight_Gsf_L1EGMT", 1); tree_->SetBranchAddress("HLT_Ele35_WPTight_Gsf_L1EGMT", &HLT_Ele35_WPTight_Gsf_L1EGMT_);
      tree_->SetBranchStatus("HLT_Ele38_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele40_WPTight_Gsf", 1); tree_->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf_);
      tree_->SetBranchStatus("HLT_Ele32_WPTight_Gsf_L1DoubleEG", 1); tree_->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_);
      tree_->SetBranchStatus("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", 1); tree_->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_);
      tree_->SetBranchStatus("HLT_HT450_Beamspot", 1); tree_->SetBranchAddress("HLT_HT450_Beamspot", &HLT_HT450_Beamspot_);
      tree_->SetBranchStatus("HLT_HT300_Beamspot", 1); tree_->SetBranchAddress("HLT_HT300_Beamspot", &HLT_HT300_Beamspot_);
      tree_->SetBranchStatus("HLT_ZeroBias_Beamspot", 1); tree_->SetBranchAddress("HLT_ZeroBias_Beamspot", &HLT_ZeroBias_Beamspot_);
      tree_->SetBranchStatus("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1", 1); tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_);
      tree_->SetBranchStatus("HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", 1); tree_->SetBranchAddress("HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", &HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_);
      tree_->SetBranchStatus("HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", 1); tree_->SetBranchAddress("HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", &HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_);
      tree_->SetBranchStatus("HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", 1); tree_->SetBranchAddress("HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1", &HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_);
      tree_->SetBranchStatus("HLT_IsoMu20", 1); tree_->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20_);
      tree_->SetBranchStatus("HLT_IsoMu24", 1); tree_->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24_);
      tree_->SetBranchStatus("HLT_IsoMu24_eta2p1", 1); tree_->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1_);
      tree_->SetBranchStatus("HLT_IsoMu27", 1); tree_->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27_);
      tree_->SetBranchStatus("HLT_IsoMu30", 1); tree_->SetBranchAddress("HLT_IsoMu30", &HLT_IsoMu30_);
      tree_->SetBranchStatus("HLT_UncorrectedJetE30_NoBPTX", 1); tree_->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX", &HLT_UncorrectedJetE30_NoBPTX_);
      tree_->SetBranchStatus("HLT_UncorrectedJetE30_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX3BX", &HLT_UncorrectedJetE30_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_UncorrectedJetE60_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_UncorrectedJetE60_NoBPTX3BX", &HLT_UncorrectedJetE60_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_UncorrectedJetE70_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_UncorrectedJetE70_NoBPTX3BX", &HLT_UncorrectedJetE70_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_L1SingleMu18", 1); tree_->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18_);
      tree_->SetBranchStatus("HLT_L1SingleMu25", 1); tree_->SetBranchAddress("HLT_L1SingleMu25", &HLT_L1SingleMu25_);
      tree_->SetBranchStatus("HLT_L2Mu10", 1); tree_->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10_);
      tree_->SetBranchStatus("HLT_L2Mu10_NoVertex_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_L2Mu10_NoVertex_NoBPTX", 1); tree_->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX_);
      tree_->SetBranchStatus("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", 1); tree_->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_);
      tree_->SetBranchStatus("HLT_L2Mu50", 1); tree_->SetBranchAddress("HLT_L2Mu50", &HLT_L2Mu50_);
      tree_->SetBranchStatus("HLT_L2Mu23NoVtx_2Cha", 1); tree_->SetBranchAddress("HLT_L2Mu23NoVtx_2Cha", &HLT_L2Mu23NoVtx_2Cha_);
      tree_->SetBranchStatus("HLT_L2Mu23NoVtx_2Cha_CosmicSeed", 1); tree_->SetBranchAddress("HLT_L2Mu23NoVtx_2Cha_CosmicSeed", &HLT_L2Mu23NoVtx_2Cha_CosmicSeed_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4", &HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4", &HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu50", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu50", &HLT_DoubleL2Mu50_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed", &HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched", &HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu23NoVtx_2Cha", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha", &HLT_DoubleL2Mu23NoVtx_2Cha_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched", &HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu25NoVtx_2Cha", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha", &HLT_DoubleL2Mu25NoVtx_2Cha_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched", &HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_);
      tree_->SetBranchStatus("HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4", 1); tree_->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4", &HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_);
      tree_->SetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", 1); tree_->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_);
      tree_->SetBranchStatus("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", 1); tree_->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_);
      tree_->SetBranchStatus("HLT_Mu25_TkMu0_Onia", 1); tree_->SetBranchAddress("HLT_Mu25_TkMu0_Onia", &HLT_Mu25_TkMu0_Onia_);
      tree_->SetBranchStatus("HLT_Mu30_TkMu0_Psi", 1); tree_->SetBranchAddress("HLT_Mu30_TkMu0_Psi", &HLT_Mu30_TkMu0_Psi_);
      tree_->SetBranchStatus("HLT_Mu30_TkMu0_Upsilon", 1); tree_->SetBranchAddress("HLT_Mu30_TkMu0_Upsilon", &HLT_Mu30_TkMu0_Upsilon_);
      tree_->SetBranchStatus("HLT_Mu20_TkMu0_Phi", 1); tree_->SetBranchAddress("HLT_Mu20_TkMu0_Phi", &HLT_Mu20_TkMu0_Phi_);
      tree_->SetBranchStatus("HLT_Mu25_TkMu0_Phi", 1); tree_->SetBranchAddress("HLT_Mu25_TkMu0_Phi", &HLT_Mu25_TkMu0_Phi_);
      tree_->SetBranchStatus("HLT_Mu12", 1); tree_->SetBranchAddress("HLT_Mu12", &HLT_Mu12_);
      tree_->SetBranchStatus("HLT_Mu15", 1); tree_->SetBranchAddress("HLT_Mu15", &HLT_Mu15_);
      tree_->SetBranchStatus("HLT_Mu20", 1); tree_->SetBranchAddress("HLT_Mu20", &HLT_Mu20_);
      tree_->SetBranchStatus("HLT_Mu27", 1); tree_->SetBranchAddress("HLT_Mu27", &HLT_Mu27_);
      tree_->SetBranchStatus("HLT_Mu50", 1); tree_->SetBranchAddress("HLT_Mu50", &HLT_Mu50_);
      tree_->SetBranchStatus("HLT_Mu55", 1); tree_->SetBranchAddress("HLT_Mu55", &HLT_Mu55_);
      tree_->SetBranchStatus("HLT_OldMu100", 1); tree_->SetBranchAddress("HLT_OldMu100", &HLT_OldMu100_);
      tree_->SetBranchStatus("HLT_TkMu100", 1); tree_->SetBranchAddress("HLT_TkMu100", &HLT_TkMu100_);
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
      tree_->SetBranchStatus("HLT_AK8PFJet15", 1); tree_->SetBranchAddress("HLT_AK8PFJet15", &HLT_AK8PFJet15_);
      tree_->SetBranchStatus("HLT_AK8PFJet25", 1); tree_->SetBranchAddress("HLT_AK8PFJet25", &HLT_AK8PFJet25_);
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
      tree_->SetBranchStatus("HLT_AK8PFJet550", 1); tree_->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550_);
      tree_->SetBranchStatus("HLT_PFJet15", 1); tree_->SetBranchAddress("HLT_PFJet15", &HLT_PFJet15_);
      tree_->SetBranchStatus("HLT_PFJet25", 1); tree_->SetBranchAddress("HLT_PFJet25", &HLT_PFJet25_);
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
      tree_->SetBranchStatus("HLT_PFJet550", 1); tree_->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550_);
      tree_->SetBranchStatus("HLT_PFJetFwd15", 1); tree_->SetBranchAddress("HLT_PFJetFwd15", &HLT_PFJetFwd15_);
      tree_->SetBranchStatus("HLT_PFJetFwd25", 1); tree_->SetBranchAddress("HLT_PFJetFwd25", &HLT_PFJetFwd25_);
      tree_->SetBranchStatus("HLT_PFJetFwd40", 1); tree_->SetBranchAddress("HLT_PFJetFwd40", &HLT_PFJetFwd40_);
      tree_->SetBranchStatus("HLT_PFJetFwd60", 1); tree_->SetBranchAddress("HLT_PFJetFwd60", &HLT_PFJetFwd60_);
      tree_->SetBranchStatus("HLT_PFJetFwd80", 1); tree_->SetBranchAddress("HLT_PFJetFwd80", &HLT_PFJetFwd80_);
      tree_->SetBranchStatus("HLT_PFJetFwd140", 1); tree_->SetBranchAddress("HLT_PFJetFwd140", &HLT_PFJetFwd140_);
      tree_->SetBranchStatus("HLT_PFJetFwd200", 1); tree_->SetBranchAddress("HLT_PFJetFwd200", &HLT_PFJetFwd200_);
      tree_->SetBranchStatus("HLT_PFJetFwd260", 1); tree_->SetBranchAddress("HLT_PFJetFwd260", &HLT_PFJetFwd260_);
      tree_->SetBranchStatus("HLT_PFJetFwd320", 1); tree_->SetBranchAddress("HLT_PFJetFwd320", &HLT_PFJetFwd320_);
      tree_->SetBranchStatus("HLT_PFJetFwd400", 1); tree_->SetBranchAddress("HLT_PFJetFwd400", &HLT_PFJetFwd400_);
      tree_->SetBranchStatus("HLT_PFJetFwd450", 1); tree_->SetBranchAddress("HLT_PFJetFwd450", &HLT_PFJetFwd450_);
      tree_->SetBranchStatus("HLT_PFJetFwd500", 1); tree_->SetBranchAddress("HLT_PFJetFwd500", &HLT_PFJetFwd500_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd15", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd15", &HLT_AK8PFJetFwd15_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd25", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd25", &HLT_AK8PFJetFwd25_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd40", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd40", &HLT_AK8PFJetFwd40_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd60", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd60", &HLT_AK8PFJetFwd60_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd80", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd80", &HLT_AK8PFJetFwd80_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd140", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd140", &HLT_AK8PFJetFwd140_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd200", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd200", &HLT_AK8PFJetFwd200_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd260", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd260", &HLT_AK8PFJetFwd260_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd320", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd320", &HLT_AK8PFJetFwd320_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd400", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd400", &HLT_AK8PFJetFwd400_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd450", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd450", &HLT_AK8PFJetFwd450_);
      tree_->SetBranchStatus("HLT_AK8PFJetFwd500", 1); tree_->SetBranchAddress("HLT_AK8PFJetFwd500", &HLT_AK8PFJetFwd500_);
      tree_->SetBranchStatus("HLT_PFHT180", 1); tree_->SetBranchAddress("HLT_PFHT180", &HLT_PFHT180_);
      tree_->SetBranchStatus("HLT_PFHT250", 1); tree_->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250_);
      tree_->SetBranchStatus("HLT_PFHT370", 1); tree_->SetBranchAddress("HLT_PFHT370", &HLT_PFHT370_);
      tree_->SetBranchStatus("HLT_PFHT430", 1); tree_->SetBranchAddress("HLT_PFHT430", &HLT_PFHT430_);
      tree_->SetBranchStatus("HLT_PFHT510", 1); tree_->SetBranchAddress("HLT_PFHT510", &HLT_PFHT510_);
      tree_->SetBranchStatus("HLT_PFHT590", 1); tree_->SetBranchAddress("HLT_PFHT590", &HLT_PFHT590_);
      tree_->SetBranchStatus("HLT_PFHT680", 1); tree_->SetBranchAddress("HLT_PFHT680", &HLT_PFHT680_);
      tree_->SetBranchStatus("HLT_PFHT780", 1); tree_->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780_);
      tree_->SetBranchStatus("HLT_PFHT890", 1); tree_->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890_);
      tree_->SetBranchStatus("HLT_PFHT1050", 1); tree_->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050_);
      tree_->SetBranchStatus("HLT_PFHT500_PFMET100_PFMHT100_IDTight", 1); tree_->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_);
      tree_->SetBranchStatus("HLT_PFHT500_PFMET110_PFMHT110_IDTight", 1); tree_->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_);
      tree_->SetBranchStatus("HLT_PFHT700_PFMET85_PFMHT85_IDTight", 1); tree_->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_);
      tree_->SetBranchStatus("HLT_PFHT700_PFMET95_PFMHT95_IDTight", 1); tree_->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_);
      tree_->SetBranchStatus("HLT_PFHT800_PFMET75_PFMHT75_IDTight", 1); tree_->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_);
      tree_->SetBranchStatus("HLT_PFHT800_PFMET85_PFMHT85_IDTight", 1); tree_->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_);
      tree_->SetBranchStatus("HLT_PFMET110_PFMHT110_IDTight", 1); tree_->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight_);
      tree_->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight", 1); tree_->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight_);
      tree_->SetBranchStatus("HLT_PFMET130_PFMHT130_IDTight", 1); tree_->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight_);
      tree_->SetBranchStatus("HLT_PFMET140_PFMHT140_IDTight", 1); tree_->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight_);
      tree_->SetBranchStatus("HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1", 1); tree_->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_);
      tree_->SetBranchStatus("HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1", 1); tree_->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_);
      tree_->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1", 1); tree_->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_);
      tree_->SetBranchStatus("HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1", 1); tree_->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_);
      tree_->SetBranchStatus("HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1", 1); tree_->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1", &HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_);
      tree_->SetBranchStatus("HLT_PFMET120_PFMHT120_IDTight_PFHT60", 1); tree_->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_);
      tree_->SetBranchStatus("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", 1); tree_->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_);
      tree_->SetBranchStatus("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", 1); tree_->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", &HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_);
      tree_->SetBranchStatus("HLT_PFMETTypeOne110_PFMHT110_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETTypeOne110_PFMHT110_IDTight", &HLT_PFMETTypeOne110_PFMHT110_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETTypeOne120_PFMHT120_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight", &HLT_PFMETTypeOne120_PFMHT120_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETTypeOne130_PFMHT130_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETTypeOne130_PFMHT130_IDTight", &HLT_PFMETTypeOne130_PFMHT130_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETTypeOne140_PFMHT140_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_);
      tree_->SetBranchStatus("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", 1); tree_->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_);
      tree_->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", 1); tree_->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_);
      tree_->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", 1); tree_->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_);
      tree_->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight", 1); tree_->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_);
      tree_->SetBranchStatus("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight", 1); tree_->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_);
      tree_->SetBranchStatus("HLT_L1ETMHadSeeds", 1); tree_->SetBranchAddress("HLT_L1ETMHadSeeds", &HLT_L1ETMHadSeeds_);
      tree_->SetBranchStatus("HLT_CaloMHT90", 1); tree_->SetBranchAddress("HLT_CaloMHT90", &HLT_CaloMHT90_);
      tree_->SetBranchStatus("HLT_CaloMET80_NotCleaned", 1); tree_->SetBranchAddress("HLT_CaloMET80_NotCleaned", &HLT_CaloMET80_NotCleaned_);
      tree_->SetBranchStatus("HLT_CaloMET90_NotCleaned", 1); tree_->SetBranchAddress("HLT_CaloMET90_NotCleaned", &HLT_CaloMET90_NotCleaned_);
      tree_->SetBranchStatus("HLT_CaloMET100_NotCleaned", 1); tree_->SetBranchAddress("HLT_CaloMET100_NotCleaned", &HLT_CaloMET100_NotCleaned_);
      tree_->SetBranchStatus("HLT_CaloMET110_NotCleaned", 1); tree_->SetBranchAddress("HLT_CaloMET110_NotCleaned", &HLT_CaloMET110_NotCleaned_);
      tree_->SetBranchStatus("HLT_CaloMET250_NotCleaned", 1); tree_->SetBranchAddress("HLT_CaloMET250_NotCleaned", &HLT_CaloMET250_NotCleaned_);
      tree_->SetBranchStatus("HLT_CaloMET70_HBHECleaned", 1); tree_->SetBranchAddress("HLT_CaloMET70_HBHECleaned", &HLT_CaloMET70_HBHECleaned_);
      tree_->SetBranchStatus("HLT_CaloMET80_HBHECleaned", 1); tree_->SetBranchAddress("HLT_CaloMET80_HBHECleaned", &HLT_CaloMET80_HBHECleaned_);
      tree_->SetBranchStatus("HLT_CaloMET90_HBHECleaned", 1); tree_->SetBranchAddress("HLT_CaloMET90_HBHECleaned", &HLT_CaloMET90_HBHECleaned_);
      tree_->SetBranchStatus("HLT_CaloMET100_HBHECleaned", 1); tree_->SetBranchAddress("HLT_CaloMET100_HBHECleaned", &HLT_CaloMET100_HBHECleaned_);
      tree_->SetBranchStatus("HLT_CaloMET250_HBHECleaned", 1); tree_->SetBranchAddress("HLT_CaloMET250_HBHECleaned", &HLT_CaloMET250_HBHECleaned_);
      tree_->SetBranchStatus("HLT_CaloMET300_HBHECleaned", 1); tree_->SetBranchAddress("HLT_CaloMET300_HBHECleaned", &HLT_CaloMET300_HBHECleaned_);
      tree_->SetBranchStatus("HLT_CaloMET350_HBHECleaned", 1); tree_->SetBranchAddress("HLT_CaloMET350_HBHECleaned", &HLT_CaloMET350_HBHECleaned_);
      tree_->SetBranchStatus("HLT_PFMET200_NotCleaned", 1); tree_->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned_);
      tree_->SetBranchStatus("HLT_PFMET200_HBHECleaned", 1); tree_->SetBranchAddress("HLT_PFMET200_HBHECleaned", &HLT_PFMET200_HBHECleaned_);
      tree_->SetBranchStatus("HLT_PFMET250_HBHECleaned", 1); tree_->SetBranchAddress("HLT_PFMET250_HBHECleaned", &HLT_PFMET250_HBHECleaned_);
      tree_->SetBranchStatus("HLT_PFMET300_HBHECleaned", 1); tree_->SetBranchAddress("HLT_PFMET300_HBHECleaned", &HLT_PFMET300_HBHECleaned_);
      tree_->SetBranchStatus("HLT_PFMET200_HBHE_BeamHaloCleaned", 1); tree_->SetBranchAddress("HLT_PFMET200_HBHE_BeamHaloCleaned", &HLT_PFMET200_HBHE_BeamHaloCleaned_);
      tree_->SetBranchStatus("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned", 1); tree_->SetBranchAddress("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_);
      tree_->SetBranchStatus("HLT_MET105_IsoTrk50", 1); tree_->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50_);
      tree_->SetBranchStatus("HLT_MET120_IsoTrk50", 1); tree_->SetBranchAddress("HLT_MET120_IsoTrk50", &HLT_MET120_IsoTrk50_);
      tree_->SetBranchStatus("HLT_SingleJet30_Mu12_SinglePFJet40", 1); tree_->SetBranchAddress("HLT_SingleJet30_Mu12_SinglePFJet40", &HLT_SingleJet30_Mu12_SinglePFJet40_);
      tree_->SetBranchStatus("HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_DoublePFJets40_CaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_DoublePFJets40_CaloBTagDeepCSV_p71", &HLT_DoublePFJets40_CaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_DoublePFJets100_CaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_DoublePFJets100_CaloBTagDeepCSV_p71", &HLT_DoublePFJets100_CaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_DoublePFJets200_CaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_DoublePFJets200_CaloBTagDeepCSV_p71", &HLT_DoublePFJets200_CaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_DoublePFJets350_CaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_DoublePFJets350_CaloBTagDeepCSV_p71", &HLT_DoublePFJets350_CaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", 1); tree_->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
      tree_->SetBranchStatus("HLT_Photon300_NoHE", 1); tree_->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", 1); tree_->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_);
      tree_->SetBranchStatus("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", 1); tree_->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_);
      tree_->SetBranchStatus("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", 1); tree_->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_);
      tree_->SetBranchStatus("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", 1); tree_->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu17_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_Mu19_TrkIsoVVL", 1); tree_->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet20_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5", &HLT_BTagMu_AK4DiJet20_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet40_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5", &HLT_BTagMu_AK4DiJet40_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet70_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5", &HLT_BTagMu_AK4DiJet70_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet110_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5", &HLT_BTagMu_AK4DiJet110_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet170_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5", &HLT_BTagMu_AK4DiJet170_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4Jet300_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5", &HLT_BTagMu_AK4Jet300_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK8DiJet170_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5", &HLT_BTagMu_AK8DiJet170_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK8Jet170_DoubleMu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK8Jet170_DoubleMu5", &HLT_BTagMu_AK8Jet170_DoubleMu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK8Jet300_Mu5", 1); tree_->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet20_Mu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5_noalgo", &HLT_BTagMu_AK4DiJet20_Mu5_noalgo_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet40_Mu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5_noalgo", &HLT_BTagMu_AK4DiJet40_Mu5_noalgo_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet70_Mu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5_noalgo", &HLT_BTagMu_AK4DiJet70_Mu5_noalgo_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet110_Mu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5_noalgo", &HLT_BTagMu_AK4DiJet110_Mu5_noalgo_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4DiJet170_Mu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5_noalgo", &HLT_BTagMu_AK4DiJet170_Mu5_noalgo_);
      tree_->SetBranchStatus("HLT_BTagMu_AK4Jet300_Mu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5_noalgo", &HLT_BTagMu_AK4Jet300_Mu5_noalgo_);
      tree_->SetBranchStatus("HLT_BTagMu_AK8DiJet170_Mu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5_noalgo", &HLT_BTagMu_AK8DiJet170_Mu5_noalgo_);
      tree_->SetBranchStatus("HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo", &HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_);
      tree_->SetBranchStatus("HLT_BTagMu_AK8Jet300_Mu5_noalgo", 1); tree_->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5_noalgo", &HLT_BTagMu_AK8Jet300_Mu5_noalgo_);
      tree_->SetBranchStatus("HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL", &HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", 1); tree_->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_);
      tree_->SetBranchStatus("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", 1); tree_->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_);
      tree_->SetBranchStatus("HLT_Mu12_DoublePhoton20", 1); tree_->SetBranchAddress("HLT_Mu12_DoublePhoton20", &HLT_Mu12_DoublePhoton20_);
      tree_->SetBranchStatus("HLT_TriplePhoton_20_20_20_CaloIdLV2", 1); tree_->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2", &HLT_TriplePhoton_20_20_20_CaloIdLV2_);
      tree_->SetBranchStatus("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL", 1); tree_->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_);
      tree_->SetBranchStatus("HLT_TriplePhoton_30_30_10_CaloIdLV2", 1); tree_->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2", &HLT_TriplePhoton_30_30_10_CaloIdLV2_);
      tree_->SetBranchStatus("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL", 1); tree_->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_);
      tree_->SetBranchStatus("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL", 1); tree_->SetBranchAddress("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_);
      tree_->SetBranchStatus("HLT_Photon20", 1); tree_->SetBranchAddress("HLT_Photon20", &HLT_Photon20_);
      tree_->SetBranchStatus("HLT_Photon33", 1); tree_->SetBranchAddress("HLT_Photon33", &HLT_Photon33_);
      tree_->SetBranchStatus("HLT_Photon50", 1); tree_->SetBranchAddress("HLT_Photon50", &HLT_Photon50_);
      tree_->SetBranchStatus("HLT_Photon75", 1); tree_->SetBranchAddress("HLT_Photon75", &HLT_Photon75_);
      tree_->SetBranchStatus("HLT_Photon90", 1); tree_->SetBranchAddress("HLT_Photon90", &HLT_Photon90_);
      tree_->SetBranchStatus("HLT_Photon120", 1); tree_->SetBranchAddress("HLT_Photon120", &HLT_Photon120_);
      tree_->SetBranchStatus("HLT_Photon150", 1); tree_->SetBranchAddress("HLT_Photon150", &HLT_Photon150_);
      tree_->SetBranchStatus("HLT_Photon175", 1); tree_->SetBranchAddress("HLT_Photon175", &HLT_Photon175_);
      tree_->SetBranchStatus("HLT_Photon200", 1); tree_->SetBranchAddress("HLT_Photon200", &HLT_Photon200_);
      tree_->SetBranchStatus("HLT_Photon100EB_TightID_TightIso", 1); tree_->SetBranchAddress("HLT_Photon100EB_TightID_TightIso", &HLT_Photon100EB_TightID_TightIso_);
      tree_->SetBranchStatus("HLT_Photon110EB_TightID_TightIso", 1); tree_->SetBranchAddress("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso_);
      tree_->SetBranchStatus("HLT_Photon120EB_TightID_TightIso", 1); tree_->SetBranchAddress("HLT_Photon120EB_TightID_TightIso", &HLT_Photon120EB_TightID_TightIso_);
      tree_->SetBranchStatus("HLT_Photon100EBHE10", 1); tree_->SetBranchAddress("HLT_Photon100EBHE10", &HLT_Photon100EBHE10_);
      tree_->SetBranchStatus("HLT_Photon100EEHE10", 1); tree_->SetBranchAddress("HLT_Photon100EEHE10", &HLT_Photon100EEHE10_);
      tree_->SetBranchStatus("HLT_Photon100EE_TightID_TightIso", 1); tree_->SetBranchAddress("HLT_Photon100EE_TightID_TightIso", &HLT_Photon100EE_TightID_TightIso_);
      tree_->SetBranchStatus("HLT_Photon50_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3", 1); tree_->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_);
      tree_->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3", 1); tree_->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_);
      tree_->SetBranchStatus("HLT_Photon90_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon120_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon165_R9Id90_HE10_IsoM", 1); tree_->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM_);
      tree_->SetBranchStatus("HLT_Photon90_CaloIdL_PFHT700", 1); tree_->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT700", &HLT_Photon90_CaloIdL_PFHT700_);
      tree_->SetBranchStatus("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", 1); tree_->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_);
      tree_->SetBranchStatus("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", 1); tree_->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_);
      tree_->SetBranchStatus("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", 1); tree_->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_);
      tree_->SetBranchStatus("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", 1); tree_->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_);
      tree_->SetBranchStatus("HLT_Photon35_TwoProngs35", 1); tree_->SetBranchAddress("HLT_Photon35_TwoProngs35", &HLT_Photon35_TwoProngs35_);
      tree_->SetBranchStatus("HLT_IsoMu24_TwoProngs35", 1); tree_->SetBranchAddress("HLT_IsoMu24_TwoProngs35", &HLT_IsoMu24_TwoProngs35_);
      tree_->SetBranchStatus("HLT_Dimuon0_Jpsi_L1_NoOS", 1); tree_->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_NoOS", &HLT_Dimuon0_Jpsi_L1_NoOS_);
      tree_->SetBranchStatus("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", 1); tree_->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", &HLT_Dimuon0_Jpsi_NoVertexing_NoOS_);
      tree_->SetBranchStatus("HLT_Dimuon0_Jpsi", 1); tree_->SetBranchAddress("HLT_Dimuon0_Jpsi", &HLT_Dimuon0_Jpsi_);
      tree_->SetBranchStatus("HLT_Dimuon0_Jpsi_NoVertexing", 1); tree_->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing", &HLT_Dimuon0_Jpsi_NoVertexing_);
      tree_->SetBranchStatus("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", 1); tree_->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_);
      tree_->SetBranchStatus("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", 1); tree_->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_);
      tree_->SetBranchStatus("HLT_Dimuon0_Jpsi3p5_Muon2", 1); tree_->SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", &HLT_Dimuon0_Jpsi3p5_Muon2_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_4p5", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5", &HLT_Dimuon0_Upsilon_L1_4p5_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_5", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5", &HLT_Dimuon0_Upsilon_L1_5_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_4p5NoOS", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5NoOS", &HLT_Dimuon0_Upsilon_L1_4p5NoOS_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_4p5er2p0", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0", &HLT_Dimuon0_Upsilon_L1_4p5er2p0_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", &HLT_Dimuon0_Upsilon_L1_4p5er2p0M_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_NoVertexing", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_NoVertexing", &HLT_Dimuon0_Upsilon_NoVertexing_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_L1_5M", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5M", &HLT_Dimuon0_Upsilon_L1_5M_);
      tree_->SetBranchStatus("HLT_Dimuon0_LowMass_L1_0er1p5R", 1); tree_->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5R", &HLT_Dimuon0_LowMass_L1_0er1p5R_);
      tree_->SetBranchStatus("HLT_Dimuon0_LowMass_L1_0er1p5", 1); tree_->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5", &HLT_Dimuon0_LowMass_L1_0er1p5_);
      tree_->SetBranchStatus("HLT_Dimuon0_LowMass", 1); tree_->SetBranchAddress("HLT_Dimuon0_LowMass", &HLT_Dimuon0_LowMass_);
      tree_->SetBranchStatus("HLT_Dimuon0_LowMass_L1_4", 1); tree_->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4", &HLT_Dimuon0_LowMass_L1_4_);
      tree_->SetBranchStatus("HLT_Dimuon0_LowMass_L1_4R", 1); tree_->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4R", &HLT_Dimuon0_LowMass_L1_4R_);
      tree_->SetBranchStatus("HLT_Dimuon0_LowMass_L1_TM530", 1); tree_->SetBranchAddress("HLT_Dimuon0_LowMass_L1_TM530", &HLT_Dimuon0_LowMass_L1_TM530_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_Muon_L1_TM0", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_L1_TM0", &HLT_Dimuon0_Upsilon_Muon_L1_TM0_);
      tree_->SetBranchStatus("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", 1); tree_->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", &HLT_Dimuon0_Upsilon_Muon_NoL1Mass_);
      tree_->SetBranchStatus("HLT_TripleMu_5_3_3_Mass3p8_DZ", 1); tree_->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8_DZ", &HLT_TripleMu_5_3_3_Mass3p8_DZ_);
      tree_->SetBranchStatus("HLT_TripleMu_10_5_5_DZ", 1); tree_->SetBranchAddress("HLT_TripleMu_10_5_5_DZ", &HLT_TripleMu_10_5_5_DZ_);
      tree_->SetBranchStatus("HLT_TripleMu_12_10_5", 1); tree_->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5_);
      tree_->SetBranchStatus("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", 1); tree_->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_);
      tree_->SetBranchStatus("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", 1); tree_->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_);
      tree_->SetBranchStatus("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", 1); tree_->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_);
      tree_->SetBranchStatus("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", 1); tree_->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_);
      tree_->SetBranchStatus("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", 1); tree_->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", &HLT_DoubleMu3_DZ_PFMET50_PFMHT60_);
      tree_->SetBranchStatus("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", 1); tree_->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", &HLT_DoubleMu3_DZ_PFMET70_PFMHT70_);
      tree_->SetBranchStatus("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", 1); tree_->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", &HLT_DoubleMu3_DZ_PFMET90_PFMHT90_);
      tree_->SetBranchStatus("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", 1); tree_->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", &HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_);
      tree_->SetBranchStatus("HLT_DoubleMu4_Jpsi_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_Jpsi_Displaced", &HLT_DoubleMu4_Jpsi_Displaced_);
      tree_->SetBranchStatus("HLT_DoubleMu4_Jpsi_NoVertexing", 1); tree_->SetBranchAddress("HLT_DoubleMu4_Jpsi_NoVertexing", &HLT_DoubleMu4_Jpsi_NoVertexing_);
      tree_->SetBranchStatus("HLT_DoubleMu4_JpsiTrkTrk_Displaced", 1); tree_->SetBranchAddress("HLT_DoubleMu4_JpsiTrkTrk_Displaced", &HLT_DoubleMu4_JpsiTrkTrk_Displaced_);
      tree_->SetBranchStatus("HLT_DoubleMu43NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_DoubleMu43NoFiltersNoVtx", &HLT_DoubleMu43NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_DoubleMu48NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_DoubleMu48NoFiltersNoVtx", &HLT_DoubleMu48NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_);
      tree_->SetBranchStatus("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", &HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_);
      tree_->SetBranchStatus("HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL", &HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_);
      tree_->SetBranchStatus("HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL", 1); tree_->SetBranchAddress("HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_);
      tree_->SetBranchStatus("HLT_DoubleMu33NoFiltersNoVtxDisplaced", 1); tree_->SetBranchAddress("HLT_DoubleMu33NoFiltersNoVtxDisplaced", &HLT_DoubleMu33NoFiltersNoVtxDisplaced_);
      tree_->SetBranchStatus("HLT_DoubleMu40NoFiltersNoVtxDisplaced", 1); tree_->SetBranchAddress("HLT_DoubleMu40NoFiltersNoVtxDisplaced", &HLT_DoubleMu40NoFiltersNoVtxDisplaced_);
      tree_->SetBranchStatus("HLT_DoubleMu20_7_Mass0to30_L1_DM4", 1); tree_->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4", &HLT_DoubleMu20_7_Mass0to30_L1_DM4_);
      tree_->SetBranchStatus("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG", 1); tree_->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG", &HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_);
      tree_->SetBranchStatus("HLT_HT425", 1); tree_->SetBranchAddress("HLT_HT425", &HLT_HT425_);
      tree_->SetBranchStatus("HLT_HT430_DisplacedDijet40_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack", &HLT_HT430_DisplacedDijet40_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT500_DisplacedDijet40_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT500_DisplacedDijet40_DisplacedTrack", &HLT_HT500_DisplacedDijet40_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT430_DisplacedDijet60_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack", &HLT_HT430_DisplacedDijet60_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT400_DisplacedDijet40_DisplacedTrack", 1); tree_->SetBranchAddress("HLT_HT400_DisplacedDijet40_DisplacedTrack", &HLT_HT400_DisplacedDijet40_DisplacedTrack_);
      tree_->SetBranchStatus("HLT_HT650_DisplacedDijet60_Inclusive", 1); tree_->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive", &HLT_HT650_DisplacedDijet60_Inclusive_);
      tree_->SetBranchStatus("HLT_HT550_DisplacedDijet60_Inclusive", 1); tree_->SetBranchAddress("HLT_HT550_DisplacedDijet60_Inclusive", &HLT_HT550_DisplacedDijet60_Inclusive_);
      tree_->SetBranchStatus("HLT_DiJet110_35_Mjj650_PFMET110", 1); tree_->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET110", &HLT_DiJet110_35_Mjj650_PFMET110_);
      tree_->SetBranchStatus("HLT_DiJet110_35_Mjj650_PFMET120", 1); tree_->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET120", &HLT_DiJet110_35_Mjj650_PFMET120_);
      tree_->SetBranchStatus("HLT_DiJet110_35_Mjj650_PFMET130", 1); tree_->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET130", &HLT_DiJet110_35_Mjj650_PFMET130_);
      tree_->SetBranchStatus("HLT_TripleJet110_35_35_Mjj650_PFMET110", 1); tree_->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET110", &HLT_TripleJet110_35_35_Mjj650_PFMET110_);
      tree_->SetBranchStatus("HLT_TripleJet110_35_35_Mjj650_PFMET120", 1); tree_->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET120", &HLT_TripleJet110_35_35_Mjj650_PFMET120_);
      tree_->SetBranchStatus("HLT_TripleJet110_35_35_Mjj650_PFMET130", 1); tree_->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET130", &HLT_TripleJet110_35_35_Mjj650_PFMET130_);
      tree_->SetBranchStatus("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", 1); tree_->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_);
      tree_->SetBranchStatus("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", 1); tree_->SetBranchAddress("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", &HLT_Ele28_eta2p1_WPTight_Gsf_HT150_);
      tree_->SetBranchStatus("HLT_Ele28_HighEta_SC20_Mass55", 1); tree_->SetBranchAddress("HLT_Ele28_HighEta_SC20_Mass55", &HLT_Ele28_HighEta_SC20_Mass55_);
      tree_->SetBranchStatus("HLT_DoubleMu20_7_Mass0to30_Photon23", 1); tree_->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_Photon23", &HLT_DoubleMu20_7_Mass0to30_Photon23_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", &HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", &HLT_Ele15_IsoVVVL_PFHT450_PFMET50_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT450", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450", &HLT_Ele15_IsoVVVL_PFHT450_);
      tree_->SetBranchStatus("HLT_Ele50_IsoVVVL_PFHT450", 1); tree_->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT450", &HLT_Ele50_IsoVVVL_PFHT450_);
      tree_->SetBranchStatus("HLT_Ele15_IsoVVVL_PFHT600", 1); tree_->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600_);
      tree_->SetBranchStatus("HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", 1); tree_->SetBranchAddress("HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_);
      tree_->SetBranchStatus("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", 1); tree_->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_);
      tree_->SetBranchStatus("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", 1); tree_->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", &HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", &HLT_Mu15_IsoVVVL_PFHT450_PFMET50_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT450", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450", &HLT_Mu15_IsoVVVL_PFHT450_);
      tree_->SetBranchStatus("HLT_Mu50_IsoVVVL_PFHT450", 1); tree_->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT450", &HLT_Mu50_IsoVVVL_PFHT450_);
      tree_->SetBranchStatus("HLT_Mu15_IsoVVVL_PFHT600", 1); tree_->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600_);
      tree_->SetBranchStatus("HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight", 1); tree_->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_);
      tree_->SetBranchStatus("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight", 1); tree_->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_);
      tree_->SetBranchStatus("HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight", 1); tree_->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_);
      tree_->SetBranchStatus("HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight", 1); tree_->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_);
      tree_->SetBranchStatus("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight", 1); tree_->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_);
      tree_->SetBranchStatus("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight", 1); tree_->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_);
      tree_->SetBranchStatus("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight", 1); tree_->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_);
      tree_->SetBranchStatus("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight", 1); tree_->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_);
      tree_->SetBranchStatus("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", 1); tree_->SetBranchAddress("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", &HLT_Dimuon10_PsiPrime_Barrel_Seagulls_);
      tree_->SetBranchStatus("HLT_Dimuon20_Jpsi_Barrel_Seagulls", 1); tree_->SetBranchAddress("HLT_Dimuon20_Jpsi_Barrel_Seagulls", &HLT_Dimuon20_Jpsi_Barrel_Seagulls_);
      tree_->SetBranchStatus("HLT_Dimuon12_Upsilon_y1p4", 1); tree_->SetBranchAddress("HLT_Dimuon12_Upsilon_y1p4", &HLT_Dimuon12_Upsilon_y1p4_);
      tree_->SetBranchStatus("HLT_Dimuon14_Phi_Barrel_Seagulls", 1); tree_->SetBranchAddress("HLT_Dimuon14_Phi_Barrel_Seagulls", &HLT_Dimuon14_Phi_Barrel_Seagulls_);
      tree_->SetBranchStatus("HLT_Dimuon18_PsiPrime", 1); tree_->SetBranchAddress("HLT_Dimuon18_PsiPrime", &HLT_Dimuon18_PsiPrime_);
      tree_->SetBranchStatus("HLT_Dimuon25_Jpsi", 1); tree_->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi_);
      tree_->SetBranchStatus("HLT_Dimuon18_PsiPrime_noCorrL1", 1); tree_->SetBranchAddress("HLT_Dimuon18_PsiPrime_noCorrL1", &HLT_Dimuon18_PsiPrime_noCorrL1_);
      tree_->SetBranchStatus("HLT_Dimuon24_Upsilon_noCorrL1", 1); tree_->SetBranchAddress("HLT_Dimuon24_Upsilon_noCorrL1", &HLT_Dimuon24_Upsilon_noCorrL1_);
      tree_->SetBranchStatus("HLT_Dimuon24_Phi_noCorrL1", 1); tree_->SetBranchAddress("HLT_Dimuon24_Phi_noCorrL1", &HLT_Dimuon24_Phi_noCorrL1_);
      tree_->SetBranchStatus("HLT_Dimuon25_Jpsi_noCorrL1", 1); tree_->SetBranchAddress("HLT_Dimuon25_Jpsi_noCorrL1", &HLT_Dimuon25_Jpsi_noCorrL1_);
      tree_->SetBranchStatus("HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8", 1); tree_->SetBranchAddress("HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8", &HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_);
      tree_->SetBranchStatus("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", 1); tree_->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_);
      tree_->SetBranchStatus("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", 1); tree_->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_);
      tree_->SetBranchStatus("HLT_DoubleIsoMu20_eta2p1", 1); tree_->SetBranchAddress("HLT_DoubleIsoMu20_eta2p1", &HLT_DoubleIsoMu20_eta2p1_);
      tree_->SetBranchStatus("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx", &HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_Mu8", 1); tree_->SetBranchAddress("HLT_Mu8", &HLT_Mu8_);
      tree_->SetBranchStatus("HLT_Mu17", 1); tree_->SetBranchAddress("HLT_Mu17", &HLT_Mu17_);
      tree_->SetBranchStatus("HLT_Mu19", 1); tree_->SetBranchAddress("HLT_Mu19", &HLT_Mu19_);
      tree_->SetBranchStatus("HLT_Mu17_Photon30_IsoCaloId", 1); tree_->SetBranchAddress("HLT_Mu17_Photon30_IsoCaloId", &HLT_Mu17_Photon30_IsoCaloId_);
      tree_->SetBranchStatus("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", 1); tree_->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30_);
      tree_->SetBranchStatus("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", 1); tree_->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_);
      tree_->SetBranchStatus("HLT_Ele115_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_Ele135_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele135_CaloIdVT_GsfTrkIdT", &HLT_Ele135_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_Ele145_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele145_CaloIdVT_GsfTrkIdT", &HLT_Ele145_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_Ele200_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele200_CaloIdVT_GsfTrkIdT", &HLT_Ele200_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_Ele250_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele250_CaloIdVT_GsfTrkIdT", &HLT_Ele250_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_Ele300_CaloIdVT_GsfTrkIdT", 1); tree_->SetBranchAddress("HLT_Ele300_CaloIdVT_GsfTrkIdT", &HLT_Ele300_CaloIdVT_GsfTrkIdT_);
      tree_->SetBranchStatus("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5", 1); tree_->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_);
      tree_->SetBranchStatus("HLT_PFHT330PT30_QuadPFJet_75_60_45_40", 1); tree_->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40_);
      tree_->SetBranchStatus("HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94", 1); tree_->SetBranchAddress("HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94", &HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_);
      tree_->SetBranchStatus("HLT_PFHT400_SixPFJet32", 1); tree_->SetBranchAddress("HLT_PFHT400_SixPFJet32", &HLT_PFHT400_SixPFJet32_);
      tree_->SetBranchStatus("HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59", 1); tree_->SetBranchAddress("HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59", &HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_);
      tree_->SetBranchStatus("HLT_PFHT450_SixPFJet36", 1); tree_->SetBranchAddress("HLT_PFHT450_SixPFJet36", &HLT_PFHT450_SixPFJet36_);
      tree_->SetBranchStatus("HLT_PFHT350", 1); tree_->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350_);
      tree_->SetBranchStatus("HLT_PFHT350MinPFJet15", 1); tree_->SetBranchAddress("HLT_PFHT350MinPFJet15", &HLT_PFHT350MinPFJet15_);
      tree_->SetBranchStatus("HLT_Photon60_R9Id90_CaloIdL_IsoL", 1); tree_->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL", &HLT_Photon60_R9Id90_CaloIdL_IsoL_);
      tree_->SetBranchStatus("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL", 1); tree_->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_);
      tree_->SetBranchStatus("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", 1); tree_->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_);
      tree_->SetBranchStatus("HLT_ECALHT800", 1); tree_->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800_);
      tree_->SetBranchStatus("HLT_DiSC30_18_EIso_AND_HE_Mass70", 1); tree_->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70_);
      tree_->SetBranchStatus("HLT_Physics", 1); tree_->SetBranchAddress("HLT_Physics", &HLT_Physics_);
      tree_->SetBranchStatus("HLT_Physics_part0", 1); tree_->SetBranchAddress("HLT_Physics_part0", &HLT_Physics_part0_);
      tree_->SetBranchStatus("HLT_Physics_part1", 1); tree_->SetBranchAddress("HLT_Physics_part1", &HLT_Physics_part1_);
      tree_->SetBranchStatus("HLT_Physics_part2", 1); tree_->SetBranchAddress("HLT_Physics_part2", &HLT_Physics_part2_);
      tree_->SetBranchStatus("HLT_Physics_part3", 1); tree_->SetBranchAddress("HLT_Physics_part3", &HLT_Physics_part3_);
      tree_->SetBranchStatus("HLT_Physics_part4", 1); tree_->SetBranchAddress("HLT_Physics_part4", &HLT_Physics_part4_);
      tree_->SetBranchStatus("HLT_Physics_part5", 1); tree_->SetBranchAddress("HLT_Physics_part5", &HLT_Physics_part5_);
      tree_->SetBranchStatus("HLT_Physics_part6", 1); tree_->SetBranchAddress("HLT_Physics_part6", &HLT_Physics_part6_);
      tree_->SetBranchStatus("HLT_Physics_part7", 1); tree_->SetBranchAddress("HLT_Physics_part7", &HLT_Physics_part7_);
      tree_->SetBranchStatus("HLT_Random", 1); tree_->SetBranchAddress("HLT_Random", &HLT_Random_);
      tree_->SetBranchStatus("HLT_ZeroBias", 1); tree_->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias_);
      tree_->SetBranchStatus("HLT_ZeroBias_Alignment", 1); tree_->SetBranchAddress("HLT_ZeroBias_Alignment", &HLT_ZeroBias_Alignment_);
      tree_->SetBranchStatus("HLT_ZeroBias_part0", 1); tree_->SetBranchAddress("HLT_ZeroBias_part0", &HLT_ZeroBias_part0_);
      tree_->SetBranchStatus("HLT_ZeroBias_part1", 1); tree_->SetBranchAddress("HLT_ZeroBias_part1", &HLT_ZeroBias_part1_);
      tree_->SetBranchStatus("HLT_ZeroBias_part2", 1); tree_->SetBranchAddress("HLT_ZeroBias_part2", &HLT_ZeroBias_part2_);
      tree_->SetBranchStatus("HLT_ZeroBias_part3", 1); tree_->SetBranchAddress("HLT_ZeroBias_part3", &HLT_ZeroBias_part3_);
      tree_->SetBranchStatus("HLT_ZeroBias_part4", 1); tree_->SetBranchAddress("HLT_ZeroBias_part4", &HLT_ZeroBias_part4_);
      tree_->SetBranchStatus("HLT_ZeroBias_part5", 1); tree_->SetBranchAddress("HLT_ZeroBias_part5", &HLT_ZeroBias_part5_);
      tree_->SetBranchStatus("HLT_ZeroBias_part6", 1); tree_->SetBranchAddress("HLT_ZeroBias_part6", &HLT_ZeroBias_part6_);
      tree_->SetBranchStatus("HLT_ZeroBias_part7", 1); tree_->SetBranchAddress("HLT_ZeroBias_part7", &HLT_ZeroBias_part7_);
      tree_->SetBranchStatus("HLT_AK4CaloJet30", 1); tree_->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30_);
      tree_->SetBranchStatus("HLT_AK4CaloJet40", 1); tree_->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40_);
      tree_->SetBranchStatus("HLT_AK4CaloJet50", 1); tree_->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50_);
      tree_->SetBranchStatus("HLT_AK4CaloJet80", 1); tree_->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80_);
      tree_->SetBranchStatus("HLT_AK4CaloJet100", 1); tree_->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100_);
      tree_->SetBranchStatus("HLT_AK4CaloJet120", 1); tree_->SetBranchAddress("HLT_AK4CaloJet120", &HLT_AK4CaloJet120_);
      tree_->SetBranchStatus("HLT_AK4PFJet30", 1); tree_->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30_);
      tree_->SetBranchStatus("HLT_AK4PFJet50", 1); tree_->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50_);
      tree_->SetBranchStatus("HLT_AK4PFJet80", 1); tree_->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80_);
      tree_->SetBranchStatus("HLT_AK4PFJet100", 1); tree_->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100_);
      tree_->SetBranchStatus("HLT_AK4PFJet120", 1); tree_->SetBranchAddress("HLT_AK4PFJet120", &HLT_AK4PFJet120_);
      tree_->SetBranchStatus("HLT_SinglePhoton10_Eta3p1ForPPRef", 1); tree_->SetBranchAddress("HLT_SinglePhoton10_Eta3p1ForPPRef", &HLT_SinglePhoton10_Eta3p1ForPPRef_);
      tree_->SetBranchStatus("HLT_SinglePhoton20_Eta3p1ForPPRef", 1); tree_->SetBranchAddress("HLT_SinglePhoton20_Eta3p1ForPPRef", &HLT_SinglePhoton20_Eta3p1ForPPRef_);
      tree_->SetBranchStatus("HLT_SinglePhoton30_Eta3p1ForPPRef", 1); tree_->SetBranchAddress("HLT_SinglePhoton30_Eta3p1ForPPRef", &HLT_SinglePhoton30_Eta3p1ForPPRef_);
      tree_->SetBranchStatus("HLT_Photon20_HoverELoose", 1); tree_->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose_);
      tree_->SetBranchStatus("HLT_Photon30_HoverELoose", 1); tree_->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose_);
      tree_->SetBranchStatus("HLT_EcalCalibration", 1); tree_->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration_);
      tree_->SetBranchStatus("HLT_HcalCalibration", 1); tree_->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration_);
      tree_->SetBranchStatus("HLT_L1UnpairedBunchBptxMinus", 1); tree_->SetBranchAddress("HLT_L1UnpairedBunchBptxMinus", &HLT_L1UnpairedBunchBptxMinus_);
      tree_->SetBranchStatus("HLT_L1UnpairedBunchBptxPlus", 1); tree_->SetBranchAddress("HLT_L1UnpairedBunchBptxPlus", &HLT_L1UnpairedBunchBptxPlus_);
      tree_->SetBranchStatus("HLT_L1NotBptxOR", 1); tree_->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR_);
      tree_->SetBranchStatus("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", 1); tree_->SetBranchAddress("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_);
      tree_->SetBranchStatus("HLT_CDC_L2cosmic_5_er1p0", 1); tree_->SetBranchAddress("HLT_CDC_L2cosmic_5_er1p0", &HLT_CDC_L2cosmic_5_er1p0_);
      tree_->SetBranchStatus("HLT_CDC_L2cosmic_5p5_er1p0", 1); tree_->SetBranchAddress("HLT_CDC_L2cosmic_5p5_er1p0", &HLT_CDC_L2cosmic_5p5_er1p0_);
      tree_->SetBranchStatus("HLT_HcalNZS", 1); tree_->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS_);
      tree_->SetBranchStatus("HLT_HcalPhiSym", 1); tree_->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym_);
      tree_->SetBranchStatus("HLT_HcalIsolatedbunch", 1); tree_->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch_);
      tree_->SetBranchStatus("HLT_IsoTrackHB", 1); tree_->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB_);
      tree_->SetBranchStatus("HLT_IsoTrackHE", 1); tree_->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE_);
      tree_->SetBranchStatus("HLT_ZeroBias_FirstCollisionAfterAbortGap", 1); tree_->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap_);
      tree_->SetBranchStatus("HLT_ZeroBias_IsolatedBunches", 1); tree_->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches_);
      tree_->SetBranchStatus("HLT_ZeroBias_FirstCollisionInTrain", 1); tree_->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain_);
      tree_->SetBranchStatus("HLT_ZeroBias_LastCollisionInTrain", 1); tree_->SetBranchAddress("HLT_ZeroBias_LastCollisionInTrain", &HLT_ZeroBias_LastCollisionInTrain_);
      tree_->SetBranchStatus("HLT_ZeroBias_FirstBXAfterTrain", 1); tree_->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain_);
      tree_->SetBranchStatus("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", 1); tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_);
      tree_->SetBranchStatus("HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1", 1); tree_->SetBranchAddress("HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_);
      tree_->SetBranchStatus("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", 1); tree_->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_);
      tree_->SetBranchStatus("HLT_Rsq0p35", 1); tree_->SetBranchAddress("HLT_Rsq0p35", &HLT_Rsq0p35_);
      tree_->SetBranchStatus("HLT_Rsq0p40", 1); tree_->SetBranchAddress("HLT_Rsq0p40", &HLT_Rsq0p40_);
      tree_->SetBranchStatus("HLT_RsqMR300_Rsq0p09_MR200", 1); tree_->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200", &HLT_RsqMR300_Rsq0p09_MR200_);
      tree_->SetBranchStatus("HLT_RsqMR320_Rsq0p09_MR200", 1); tree_->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200", &HLT_RsqMR320_Rsq0p09_MR200_);
      tree_->SetBranchStatus("HLT_RsqMR300_Rsq0p09_MR200_4jet", 1); tree_->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200_4jet", &HLT_RsqMR300_Rsq0p09_MR200_4jet_);
      tree_->SetBranchStatus("HLT_RsqMR320_Rsq0p09_MR200_4jet", 1); tree_->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200_4jet", &HLT_RsqMR320_Rsq0p09_MR200_4jet_);
      tree_->SetBranchStatus("HLT_IsoMu27_MET90", 1); tree_->SetBranchAddress("HLT_IsoMu27_MET90", &HLT_IsoMu27_MET90_);
      tree_->SetBranchStatus("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg", 1); tree_->SetBranchAddress("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_);
      tree_->SetBranchStatus("HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1", 1); tree_->SetBranchAddress("HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1", &HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_);
      tree_->SetBranchStatus("HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1", 1); tree_->SetBranchAddress("HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1", &HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_);
      tree_->SetBranchStatus("HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1", 1); tree_->SetBranchAddress("HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1", &HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_);
      tree_->SetBranchStatus("HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50", 1); tree_->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50", &HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_);
      tree_->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", 1); tree_->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_);
      tree_->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3", 1); tree_->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_);
      tree_->SetBranchStatus("HLT_PFMET100_PFMHT100_IDTight_PFHT60", 1); tree_->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_PFHT60", &HLT_PFMET100_PFMHT100_IDTight_PFHT60_);
      tree_->SetBranchStatus("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", 1); tree_->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_);
      tree_->SetBranchStatus("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", 1); tree_->SetBranchAddress("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", &HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_);
      tree_->SetBranchStatus("HLT_Mu18_Mu9_SameSign", 1); tree_->SetBranchAddress("HLT_Mu18_Mu9_SameSign", &HLT_Mu18_Mu9_SameSign_);
      tree_->SetBranchStatus("HLT_Mu18_Mu9_SameSign_DZ", 1); tree_->SetBranchAddress("HLT_Mu18_Mu9_SameSign_DZ", &HLT_Mu18_Mu9_SameSign_DZ_);
      tree_->SetBranchStatus("HLT_Mu18_Mu9", 1); tree_->SetBranchAddress("HLT_Mu18_Mu9", &HLT_Mu18_Mu9_);
      tree_->SetBranchStatus("HLT_Mu18_Mu9_DZ", 1); tree_->SetBranchAddress("HLT_Mu18_Mu9_DZ", &HLT_Mu18_Mu9_DZ_);
      tree_->SetBranchStatus("HLT_Mu20_Mu10_SameSign", 1); tree_->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign_);
      tree_->SetBranchStatus("HLT_Mu20_Mu10_SameSign_DZ", 1); tree_->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ_);
      tree_->SetBranchStatus("HLT_Mu20_Mu10", 1); tree_->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10_);
      tree_->SetBranchStatus("HLT_Mu20_Mu10_DZ", 1); tree_->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ_);
      tree_->SetBranchStatus("HLT_Mu23_Mu12_SameSign", 1); tree_->SetBranchAddress("HLT_Mu23_Mu12_SameSign", &HLT_Mu23_Mu12_SameSign_);
      tree_->SetBranchStatus("HLT_Mu23_Mu12_SameSign_DZ", 1); tree_->SetBranchAddress("HLT_Mu23_Mu12_SameSign_DZ", &HLT_Mu23_Mu12_SameSign_DZ_);
      tree_->SetBranchStatus("HLT_Mu23_Mu12", 1); tree_->SetBranchAddress("HLT_Mu23_Mu12", &HLT_Mu23_Mu12_);
      tree_->SetBranchStatus("HLT_Mu23_Mu12_DZ", 1); tree_->SetBranchAddress("HLT_Mu23_Mu12_DZ", &HLT_Mu23_Mu12_DZ_);
      tree_->SetBranchStatus("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", 1); tree_->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_);
      tree_->SetBranchStatus("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", 1); tree_->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", &HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_);
      tree_->SetBranchStatus("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", 1); tree_->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60_);
      tree_->SetBranchStatus("HLT_TripleMu_5_3_3_Mass3p8_DCA", 1); tree_->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8_DCA", &HLT_TripleMu_5_3_3_Mass3p8_DCA_);
      tree_->SetBranchStatus("HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", 1); tree_->SetBranchAddress("HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_);
      tree_->SetBranchStatus("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", 1); tree_->SetBranchAddress("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_);
      tree_->SetBranchStatus("HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", 1); tree_->SetBranchAddress("HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_);
      tree_->SetBranchStatus("HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2", 1); tree_->SetBranchAddress("HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_);
      tree_->SetBranchStatus("HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2", 1); tree_->SetBranchAddress("HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_);
      tree_->SetBranchStatus("HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2", 1); tree_->SetBranchAddress("HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_);
      tree_->SetBranchStatus("HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2", 1); tree_->SetBranchAddress("HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2", &HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_);
      tree_->SetBranchStatus("HLT_QuadPFJet98_83_71_15", 1); tree_->SetBranchAddress("HLT_QuadPFJet98_83_71_15", &HLT_QuadPFJet98_83_71_15_);
      tree_->SetBranchStatus("HLT_QuadPFJet103_88_75_15", 1); tree_->SetBranchAddress("HLT_QuadPFJet103_88_75_15", &HLT_QuadPFJet103_88_75_15_);
      tree_->SetBranchStatus("HLT_QuadPFJet105_88_76_15", 1); tree_->SetBranchAddress("HLT_QuadPFJet105_88_76_15", &HLT_QuadPFJet105_88_76_15_);
      tree_->SetBranchStatus("HLT_QuadPFJet111_90_80_15", 1); tree_->SetBranchAddress("HLT_QuadPFJet111_90_80_15", &HLT_QuadPFJet111_90_80_15_);
      tree_->SetBranchStatus("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17", 1); tree_->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17", &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_);
      tree_->SetBranchStatus("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1", 1); tree_->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1", &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_);
      tree_->SetBranchStatus("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02", 1); tree_->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_);
      tree_->SetBranchStatus("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2", 1); tree_->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_);
      tree_->SetBranchStatus("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4", 1); tree_->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_);
      tree_->SetBranchStatus("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55", 1); tree_->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_);
      tree_->SetBranchStatus("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto", 1); tree_->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_);
      tree_->SetBranchStatus("HLT_Mu12_IP6_part0", 1); tree_->SetBranchAddress("HLT_Mu12_IP6_part0", &HLT_Mu12_IP6_part0_);
      tree_->SetBranchStatus("HLT_Mu12_IP6_part1", 1); tree_->SetBranchAddress("HLT_Mu12_IP6_part1", &HLT_Mu12_IP6_part1_);
      tree_->SetBranchStatus("HLT_Mu12_IP6_part2", 1); tree_->SetBranchAddress("HLT_Mu12_IP6_part2", &HLT_Mu12_IP6_part2_);
      tree_->SetBranchStatus("HLT_Mu12_IP6_part3", 1); tree_->SetBranchAddress("HLT_Mu12_IP6_part3", &HLT_Mu12_IP6_part3_);
      tree_->SetBranchStatus("HLT_Mu12_IP6_part4", 1); tree_->SetBranchAddress("HLT_Mu12_IP6_part4", &HLT_Mu12_IP6_part4_);
      tree_->SetBranchStatus("HLT_Mu9_IP5_part0", 1); tree_->SetBranchAddress("HLT_Mu9_IP5_part0", &HLT_Mu9_IP5_part0_);
      tree_->SetBranchStatus("HLT_Mu9_IP5_part1", 1); tree_->SetBranchAddress("HLT_Mu9_IP5_part1", &HLT_Mu9_IP5_part1_);
      tree_->SetBranchStatus("HLT_Mu9_IP5_part2", 1); tree_->SetBranchAddress("HLT_Mu9_IP5_part2", &HLT_Mu9_IP5_part2_);
      tree_->SetBranchStatus("HLT_Mu9_IP5_part3", 1); tree_->SetBranchAddress("HLT_Mu9_IP5_part3", &HLT_Mu9_IP5_part3_);
      tree_->SetBranchStatus("HLT_Mu9_IP5_part4", 1); tree_->SetBranchAddress("HLT_Mu9_IP5_part4", &HLT_Mu9_IP5_part4_);
      tree_->SetBranchStatus("HLT_Mu7_IP4_part0", 1); tree_->SetBranchAddress("HLT_Mu7_IP4_part0", &HLT_Mu7_IP4_part0_);
      tree_->SetBranchStatus("HLT_Mu7_IP4_part1", 1); tree_->SetBranchAddress("HLT_Mu7_IP4_part1", &HLT_Mu7_IP4_part1_);
      tree_->SetBranchStatus("HLT_Mu7_IP4_part2", 1); tree_->SetBranchAddress("HLT_Mu7_IP4_part2", &HLT_Mu7_IP4_part2_);
      tree_->SetBranchStatus("HLT_Mu7_IP4_part3", 1); tree_->SetBranchAddress("HLT_Mu7_IP4_part3", &HLT_Mu7_IP4_part3_);
      tree_->SetBranchStatus("HLT_Mu7_IP4_part4", 1); tree_->SetBranchAddress("HLT_Mu7_IP4_part4", &HLT_Mu7_IP4_part4_);
      tree_->SetBranchStatus("HLT_Mu9_IP4_part0", 1); tree_->SetBranchAddress("HLT_Mu9_IP4_part0", &HLT_Mu9_IP4_part0_);
      tree_->SetBranchStatus("HLT_Mu9_IP4_part1", 1); tree_->SetBranchAddress("HLT_Mu9_IP4_part1", &HLT_Mu9_IP4_part1_);
      tree_->SetBranchStatus("HLT_Mu9_IP4_part2", 1); tree_->SetBranchAddress("HLT_Mu9_IP4_part2", &HLT_Mu9_IP4_part2_);
      tree_->SetBranchStatus("HLT_Mu9_IP4_part3", 1); tree_->SetBranchAddress("HLT_Mu9_IP4_part3", &HLT_Mu9_IP4_part3_);
      tree_->SetBranchStatus("HLT_Mu9_IP4_part4", 1); tree_->SetBranchAddress("HLT_Mu9_IP4_part4", &HLT_Mu9_IP4_part4_);
      tree_->SetBranchStatus("HLT_Mu8_IP5_part0", 1); tree_->SetBranchAddress("HLT_Mu8_IP5_part0", &HLT_Mu8_IP5_part0_);
      tree_->SetBranchStatus("HLT_Mu8_IP5_part1", 1); tree_->SetBranchAddress("HLT_Mu8_IP5_part1", &HLT_Mu8_IP5_part1_);
      tree_->SetBranchStatus("HLT_Mu8_IP5_part2", 1); tree_->SetBranchAddress("HLT_Mu8_IP5_part2", &HLT_Mu8_IP5_part2_);
      tree_->SetBranchStatus("HLT_Mu8_IP5_part3", 1); tree_->SetBranchAddress("HLT_Mu8_IP5_part3", &HLT_Mu8_IP5_part3_);
      tree_->SetBranchStatus("HLT_Mu8_IP5_part4", 1); tree_->SetBranchAddress("HLT_Mu8_IP5_part4", &HLT_Mu8_IP5_part4_);
      tree_->SetBranchStatus("HLT_Mu8_IP6_part0", 1); tree_->SetBranchAddress("HLT_Mu8_IP6_part0", &HLT_Mu8_IP6_part0_);
      tree_->SetBranchStatus("HLT_Mu8_IP6_part1", 1); tree_->SetBranchAddress("HLT_Mu8_IP6_part1", &HLT_Mu8_IP6_part1_);
      tree_->SetBranchStatus("HLT_Mu8_IP6_part2", 1); tree_->SetBranchAddress("HLT_Mu8_IP6_part2", &HLT_Mu8_IP6_part2_);
      tree_->SetBranchStatus("HLT_Mu8_IP6_part3", 1); tree_->SetBranchAddress("HLT_Mu8_IP6_part3", &HLT_Mu8_IP6_part3_);
      tree_->SetBranchStatus("HLT_Mu8_IP6_part4", 1); tree_->SetBranchAddress("HLT_Mu8_IP6_part4", &HLT_Mu8_IP6_part4_);
      tree_->SetBranchStatus("HLT_Mu9_IP6_part0", 1); tree_->SetBranchAddress("HLT_Mu9_IP6_part0", &HLT_Mu9_IP6_part0_);
      tree_->SetBranchStatus("HLT_Mu9_IP6_part1", 1); tree_->SetBranchAddress("HLT_Mu9_IP6_part1", &HLT_Mu9_IP6_part1_);
      tree_->SetBranchStatus("HLT_Mu9_IP6_part2", 1); tree_->SetBranchAddress("HLT_Mu9_IP6_part2", &HLT_Mu9_IP6_part2_);
      tree_->SetBranchStatus("HLT_Mu9_IP6_part3", 1); tree_->SetBranchAddress("HLT_Mu9_IP6_part3", &HLT_Mu9_IP6_part3_);
      tree_->SetBranchStatus("HLT_Mu9_IP6_part4", 1); tree_->SetBranchAddress("HLT_Mu9_IP6_part4", &HLT_Mu9_IP6_part4_);
      tree_->SetBranchStatus("HLT_Mu8_IP3_part0", 1); tree_->SetBranchAddress("HLT_Mu8_IP3_part0", &HLT_Mu8_IP3_part0_);
      tree_->SetBranchStatus("HLT_Mu8_IP3_part1", 1); tree_->SetBranchAddress("HLT_Mu8_IP3_part1", &HLT_Mu8_IP3_part1_);
      tree_->SetBranchStatus("HLT_Mu8_IP3_part2", 1); tree_->SetBranchAddress("HLT_Mu8_IP3_part2", &HLT_Mu8_IP3_part2_);
      tree_->SetBranchStatus("HLT_Mu8_IP3_part3", 1); tree_->SetBranchAddress("HLT_Mu8_IP3_part3", &HLT_Mu8_IP3_part3_);
      tree_->SetBranchStatus("HLT_Mu8_IP3_part4", 1); tree_->SetBranchAddress("HLT_Mu8_IP3_part4", &HLT_Mu8_IP3_part4_);
      tree_->SetBranchStatus("HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", 1); tree_->SetBranchAddress("HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_);
      tree_->SetBranchStatus("HLT_TrkMu6NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_TrkMu6NoFiltersNoVtx", &HLT_TrkMu6NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_TrkMu16NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_TrkMu16NoFiltersNoVtx", &HLT_TrkMu16NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLT_DoubleTrkMu_16_6_NoFiltersNoVtx", 1); tree_->SetBranchAddress("HLT_DoubleTrkMu_16_6_NoFiltersNoVtx", &HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_);
      tree_->SetBranchStatus("HLTriggerFinalPath", 1); tree_->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath);
      are_HLT_loaded_ = true;
      tree_->GetEntry(current_entry_);
    }
  }
  
  void loadSvs(){
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
  
  void loadElectrons(){
    if(!are_Electron_loaded_){
      tree_->SetBranchStatus("Electron_deltaEtaSC", 1); tree_->SetBranchAddress("Electron_deltaEtaSC", &Electron_deltaEtaSC_);
      tree_->SetBranchStatus("Electron_dr03EcalRecHitSumEt", 1); tree_->SetBranchAddress("Electron_dr03EcalRecHitSumEt", &Electron_dr03EcalRecHitSumEt_);
      tree_->SetBranchStatus("Electron_dr03HcalDepth1TowerSumEt", 1); tree_->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", &Electron_dr03HcalDepth1TowerSumEt_);
      tree_->SetBranchStatus("Electron_dr03TkSumPt", 1); tree_->SetBranchAddress("Electron_dr03TkSumPt", &Electron_dr03TkSumPt_);
      tree_->SetBranchStatus("Electron_dxy", 1); tree_->SetBranchAddress("Electron_dxy", &Electron_dxy_);
      tree_->SetBranchStatus("Electron_dxyErr", 1); tree_->SetBranchAddress("Electron_dxyErr", &Electron_dxyErr_);
      tree_->SetBranchStatus("Electron_dz", 1); tree_->SetBranchAddress("Electron_dz", &Electron_dz_);
      tree_->SetBranchStatus("Electron_dzErr", 1); tree_->SetBranchAddress("Electron_dzErr", &Electron_dzErr_);
      tree_->SetBranchStatus("Electron_eInvMinusPInv", 1); tree_->SetBranchAddress("Electron_eInvMinusPInv", &Electron_eInvMinusPInv_);
      tree_->SetBranchStatus("Electron_energyErr", 1); tree_->SetBranchAddress("Electron_energyErr", &Electron_energyErr_);
      tree_->SetBranchStatus("Electron_eta", 1); tree_->SetBranchAddress("Electron_eta", &Electron_eta_);
      tree_->SetBranchStatus("Electron_hoe", 1); tree_->SetBranchAddress("Electron_hoe", &Electron_hoe_);
      tree_->SetBranchStatus("Electron_ip3d", 1); tree_->SetBranchAddress("Electron_ip3d", &Electron_ip3d_);
      tree_->SetBranchStatus("Electron_mass", 1); tree_->SetBranchAddress("Electron_mass", &Electron_mass_);
      tree_->SetBranchStatus("Electron_miniPFRelIso_all", 1); tree_->SetBranchAddress("Electron_miniPFRelIso_all", &Electron_miniPFRelIso_all_);
      tree_->SetBranchStatus("Electron_miniPFRelIso_chg", 1); tree_->SetBranchAddress("Electron_miniPFRelIso_chg", &Electron_miniPFRelIso_chg_);
      tree_->SetBranchStatus("Electron_mvaFall17V1Iso", 1); tree_->SetBranchAddress("Electron_mvaFall17V1Iso", &Electron_mvaFall17V1Iso_);
      tree_->SetBranchStatus("Electron_mvaFall17V1noIso", 1); tree_->SetBranchAddress("Electron_mvaFall17V1noIso", &Electron_mvaFall17V1noIso_);
      tree_->SetBranchStatus("Electron_mvaFall17V2Iso", 1); tree_->SetBranchAddress("Electron_mvaFall17V2Iso", &Electron_mvaFall17V2Iso_);
      tree_->SetBranchStatus("Electron_mvaFall17V2noIso", 1); tree_->SetBranchAddress("Electron_mvaFall17V2noIso", &Electron_mvaFall17V2noIso_);
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
      tree_->SetBranchStatus("Electron_cutBased_Fall17_V1", 1); tree_->SetBranchAddress("Electron_cutBased_Fall17_V1", &Electron_cutBased_Fall17_V1_);
      tree_->SetBranchStatus("Electron_jetIdx", 1); tree_->SetBranchAddress("Electron_jetIdx", &Electron_jetIdx_);
      tree_->SetBranchStatus("Electron_pdgId", 1); tree_->SetBranchAddress("Electron_pdgId", &Electron_pdgId_);
      tree_->SetBranchStatus("Electron_photonIdx", 1); tree_->SetBranchAddress("Electron_photonIdx", &Electron_photonIdx_);
      tree_->SetBranchStatus("Electron_tightCharge", 1); tree_->SetBranchAddress("Electron_tightCharge", &Electron_tightCharge_);
      tree_->SetBranchStatus("Electron_vidNestedWPBitmap", 1); tree_->SetBranchAddress("Electron_vidNestedWPBitmap", &Electron_vidNestedWPBitmap_);
      tree_->SetBranchStatus("Electron_convVeto", 1); tree_->SetBranchAddress("Electron_convVeto", &Electron_convVeto_);
      tree_->SetBranchStatus("Electron_cutBased_HEEP", 1); tree_->SetBranchAddress("Electron_cutBased_HEEP", &Electron_cutBased_HEEP_);
      tree_->SetBranchStatus("Electron_isPFcand", 1); tree_->SetBranchAddress("Electron_isPFcand", &Electron_isPFcand_);
      tree_->SetBranchStatus("Electron_lostHits", 1); tree_->SetBranchAddress("Electron_lostHits", &Electron_lostHits_);
      tree_->SetBranchStatus("Electron_mvaFall17V1Iso_WP80", 1); tree_->SetBranchAddress("Electron_mvaFall17V1Iso_WP80", &Electron_mvaFall17V1Iso_WP80_);
      tree_->SetBranchStatus("Electron_mvaFall17V1Iso_WP90", 1); tree_->SetBranchAddress("Electron_mvaFall17V1Iso_WP90", &Electron_mvaFall17V1Iso_WP90_);
      tree_->SetBranchStatus("Electron_mvaFall17V1Iso_WPL", 1); tree_->SetBranchAddress("Electron_mvaFall17V1Iso_WPL", &Electron_mvaFall17V1Iso_WPL_);
      tree_->SetBranchStatus("Electron_mvaFall17V1noIso_WP80", 1); tree_->SetBranchAddress("Electron_mvaFall17V1noIso_WP80", &Electron_mvaFall17V1noIso_WP80_);
      tree_->SetBranchStatus("Electron_mvaFall17V1noIso_WP90", 1); tree_->SetBranchAddress("Electron_mvaFall17V1noIso_WP90", &Electron_mvaFall17V1noIso_WP90_);
      tree_->SetBranchStatus("Electron_mvaFall17V1noIso_WPL", 1); tree_->SetBranchAddress("Electron_mvaFall17V1noIso_WPL", &Electron_mvaFall17V1noIso_WPL_);
      tree_->SetBranchStatus("Electron_mvaFall17V2Iso_WP80", 1); tree_->SetBranchAddress("Electron_mvaFall17V2Iso_WP80", &Electron_mvaFall17V2Iso_WP80_);
      tree_->SetBranchStatus("Electron_mvaFall17V2Iso_WP90", 1); tree_->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", &Electron_mvaFall17V2Iso_WP90_);
      tree_->SetBranchStatus("Electron_mvaFall17V2Iso_WPL", 1); tree_->SetBranchAddress("Electron_mvaFall17V2Iso_WPL", &Electron_mvaFall17V2Iso_WPL_);
      tree_->SetBranchStatus("Electron_mvaFall17V2noIso_WP80", 1); tree_->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", &Electron_mvaFall17V2noIso_WP80_);
      tree_->SetBranchStatus("Electron_mvaFall17V2noIso_WP90", 1); tree_->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", &Electron_mvaFall17V2noIso_WP90_);
      tree_->SetBranchStatus("Electron_mvaFall17V2noIso_WPL", 1); tree_->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", &Electron_mvaFall17V2noIso_WPL_);
      tree_->SetBranchStatus("Electron_genPartIdx", 1); tree_->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx_);
      tree_->SetBranchStatus("Electron_genPartFlav", 1); tree_->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav_);
      tree_->SetBranchStatus("Electron_cleanmask", 1); tree_->SetBranchAddress("Electron_cleanmask", &Electron_cleanmask_);
      are_Electron_loaded_ = true;
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
  
  void loadFatjets(){
    if(!are_FatJet_loaded_){
      tree_->SetBranchStatus("FatJet_area", 1); tree_->SetBranchAddress("FatJet_area", &FatJet_area_);
      tree_->SetBranchStatus("FatJet_btagCMVA", 1); tree_->SetBranchAddress("FatJet_btagCMVA", &FatJet_btagCMVA_);
      tree_->SetBranchStatus("FatJet_btagCSVV2", 1); tree_->SetBranchAddress("FatJet_btagCSVV2", &FatJet_btagCSVV2_);
      tree_->SetBranchStatus("FatJet_btagDeepB", 1); tree_->SetBranchAddress("FatJet_btagDeepB", &FatJet_btagDeepB_);
      tree_->SetBranchStatus("FatJet_btagHbb", 1); tree_->SetBranchAddress("FatJet_btagHbb", &FatJet_btagHbb_);
      tree_->SetBranchStatus("FatJet_eta", 1); tree_->SetBranchAddress("FatJet_eta", &FatJet_eta_);
      tree_->SetBranchStatus("FatJet_mass", 1); tree_->SetBranchAddress("FatJet_mass", &FatJet_mass_);
      tree_->SetBranchStatus("FatJet_msoftdrop", 1); tree_->SetBranchAddress("FatJet_msoftdrop", &FatJet_msoftdrop_);
      tree_->SetBranchStatus("FatJet_n2b1", 1); tree_->SetBranchAddress("FatJet_n2b1", &FatJet_n2b1_);
      tree_->SetBranchStatus("FatJet_n3b1", 1); tree_->SetBranchAddress("FatJet_n3b1", &FatJet_n3b1_);
      tree_->SetBranchStatus("FatJet_phi", 1); tree_->SetBranchAddress("FatJet_phi", &FatJet_phi_);
      tree_->SetBranchStatus("FatJet_pt", 1); tree_->SetBranchAddress("FatJet_pt", &FatJet_pt_);
      tree_->SetBranchStatus("FatJet_rawFactor", 1); tree_->SetBranchAddress("FatJet_rawFactor", &FatJet_rawFactor_);
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
  
  void loadSubjets(){
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
  
  void loadSubgenjetak8s(){
    if(!are_SubGenJetAK8_loaded_){
      tree_->SetBranchStatus("SubGenJetAK8_eta", 1); tree_->SetBranchAddress("SubGenJetAK8_eta", &SubGenJetAK8_eta_);
      tree_->SetBranchStatus("SubGenJetAK8_mass", 1); tree_->SetBranchAddress("SubGenJetAK8_mass", &SubGenJetAK8_mass_);
      tree_->SetBranchStatus("SubGenJetAK8_phi", 1); tree_->SetBranchAddress("SubGenJetAK8_phi", &SubGenJetAK8_phi_);
      tree_->SetBranchStatus("SubGenJetAK8_pt", 1); tree_->SetBranchAddress("SubGenJetAK8_pt", &SubGenJetAK8_pt_);
      are_SubGenJetAK8_loaded_ = true;
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
  

  const vector<Jets>& jets(){
    if( Jet_.size() > 0) return Jet_;
    loadJets();
  	Jet_.reserve(nJet);
    for(size_t idx = 0; idx < nJet; ++idx ){
      Jets obj;
      obj.setarea(Jet_area_[idx]);
      obj.setbtagCMVA(Jet_btagCMVA_[idx]);
      obj.setbtagCSVV2(Jet_btagCSVV2_[idx]);
      obj.setbtagDeepB(Jet_btagDeepB_[idx]);
      obj.setbtagDeepC(Jet_btagDeepC_[idx]);
      obj.setbtagDeepFlavB(Jet_btagDeepFlavB_[idx]);
      obj.setchEmEF(Jet_chEmEF_[idx]);
      obj.setchHEF(Jet_chHEF_[idx]);
      obj.setmuEF(Jet_muEF_[idx]);
      obj.setneEmEF(Jet_neEmEF_[idx]);
      obj.setneHEF(Jet_neHEF_[idx]);
      obj.setqgl(Jet_qgl_[idx]);
      obj.setrawFactor(Jet_rawFactor_[idx]);
      obj.setbRegCorr(Jet_bRegCorr_[idx]);
      obj.setbRegRes(Jet_bRegRes_[idx]);
      obj.setelectronIdx1(Jet_electronIdx1_[idx]);
      obj.setelectronIdx2(Jet_electronIdx2_[idx]);
      obj.setjetId(Jet_jetId_[idx]);
      obj.setmuonIdx1(Jet_muonIdx1_[idx]);
      obj.setmuonIdx2(Jet_muonIdx2_[idx]);
      obj.setnConstituents(Jet_nConstituents_[idx]);
      obj.setnElectrons(Jet_nElectrons_[idx]);
      obj.setnMuons(Jet_nMuons_[idx]);
      obj.setpuId(Jet_puId_[idx]);
      obj.setgenJetIdx(Jet_genJetIdx_[idx]);
      obj.sethadronFlavour(Jet_hadronFlavour_[idx]);
      obj.setpartonFlavour(Jet_partonFlavour_[idx]);
      obj.setcleanmask(Jet_cleanmask_[idx]);
      obj.setLorentzVector(Jet_pt_[idx], Jet_eta_[idx], Jet_phi_[idx], Jet_mass_[idx]);
      Jet_.push_back( obj );
    }
    return Jet_;
  }
  
  const vector<Isotracks>& isotracks(){
    if( IsoTrack_.size() > 0) return IsoTrack_;
    loadIsotracks();
  	IsoTrack_.reserve(nIsoTrack);
    for(size_t idx = 0; idx < nIsoTrack; ++idx ){
      Isotracks obj;
      obj.setdxy(IsoTrack_dxy_[idx]);
      obj.setdz(IsoTrack_dz_[idx]);
      obj.setpfRelIso03_all(IsoTrack_pfRelIso03_all_[idx]);
      obj.setpfRelIso03_chg(IsoTrack_pfRelIso03_chg_[idx]);
      obj.setminiPFRelIso_all(IsoTrack_miniPFRelIso_all_[idx]);
      obj.setminiPFRelIso_chg(IsoTrack_miniPFRelIso_chg_[idx]);
      obj.setpdgId(IsoTrack_pdgId_[idx]);
      obj.setisHighPurityTrack(IsoTrack_isHighPurityTrack_[idx]);
      obj.setisPFcand(IsoTrack_isPFcand_[idx]);
      obj.setLorentzVector(IsoTrack_pt_[idx], IsoTrack_eta_[idx], IsoTrack_phi_[idx]);
      IsoTrack_.push_back( obj );
    }
    return IsoTrack_;
  }
  
  const vector<Genjetak8s>& genjetak8s(){
    if( GenJetAK8_.size() > 0) return GenJetAK8_;
    loadGenjetak8s();
  	GenJetAK8_.reserve(nGenJetAK8);
    for(size_t idx = 0; idx < nGenJetAK8; ++idx ){
      Genjetak8s obj;
      obj.setpartonFlavour(GenJetAK8_partonFlavour_[idx]);
      obj.sethadronFlavour(GenJetAK8_hadronFlavour_[idx]);
      obj.setLorentzVector(GenJetAK8_pt_[idx], GenJetAK8_eta_[idx], GenJetAK8_phi_[idx], GenJetAK8_mass_[idx]);
      GenJetAK8_.push_back( obj );
    }
    return GenJetAK8_;
  }
  
  const vector<Genvistaus>& genvistaus(){
    if( GenVisTau_.size() > 0) return GenVisTau_;
    loadGenvistaus();
  	GenVisTau_.reserve(nGenVisTau);
    for(size_t idx = 0; idx < nGenVisTau; ++idx ){
      Genvistaus obj;
      obj.setcharge(GenVisTau_charge_[idx]);
      obj.setgenPartIdxMother(GenVisTau_genPartIdxMother_[idx]);
      obj.setstatus(GenVisTau_status_[idx]);
      obj.setLorentzVector(GenVisTau_pt_[idx], GenVisTau_eta_[idx], GenVisTau_phi_[idx], GenVisTau_mass_[idx]);
      GenVisTau_.push_back( obj );
    }
    return GenVisTau_;
  }
  
  const Calomet calomet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadCalomet();
  
    Calomet obj;
    obj.setsumEt(CaloMET_sumEt_);
    
  
    return obj;
  }
  
  const vector<Gendressedleptons>& gendressedleptons(){
    if( GenDressedLepton_.size() > 0) return GenDressedLepton_;
    loadGendressedleptons();
  	GenDressedLepton_.reserve(nGenDressedLepton);
    for(size_t idx = 0; idx < nGenDressedLepton; ++idx ){
      Gendressedleptons obj;
      obj.setpdgId(GenDressedLepton_pdgId_[idx]);
      obj.setLorentzVector(GenDressedLepton_pt_[idx], GenDressedLepton_eta_[idx], GenDressedLepton_phi_[idx], GenDressedLepton_mass_[idx]);
      GenDressedLepton_.push_back( obj );
    }
    return GenDressedLepton_;
  }
  
  const vector<Lheparts>& lheparts(){
    if( LHEPart_.size() > 0) return LHEPart_;
    loadLheparts();
  	LHEPart_.reserve(nLHEPart);
    for(size_t idx = 0; idx < nLHEPart; ++idx ){
      Lheparts obj;
      obj.setpdgId(LHEPart_pdgId_[idx]);
      obj.setLorentzVector(LHEPart_pt_[idx], LHEPart_eta_[idx], LHEPart_phi_[idx], LHEPart_mass_[idx]);
      LHEPart_.push_back( obj );
    }
    return LHEPart_;
  }
  
  const Pv pv(){
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
  
  const Generator generator(){
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
  
  const vector<Trigobjs>& trigobjs(){
    if( TrigObj_.size() > 0) return TrigObj_;
    loadTrigobjs();
  	TrigObj_.reserve(nTrigObj);
    for(size_t idx = 0; idx < nTrigObj; ++idx ){
      Trigobjs obj;
      obj.setl1pt(TrigObj_l1pt_[idx]);
      obj.setl1pt_2(TrigObj_l1pt_2_[idx]);
      obj.setl2pt(TrigObj_l2pt_[idx]);
      obj.setid(TrigObj_id_[idx]);
      obj.setl1iso(TrigObj_l1iso_[idx]);
      obj.setl1charge(TrigObj_l1charge_[idx]);
      obj.setfilterBits(TrigObj_filterBits_[idx]);
      obj.setLorentzVector(TrigObj_pt_[idx], TrigObj_eta_[idx], TrigObj_phi_[idx]);
      TrigObj_.push_back( obj );
    }
    return TrigObj_;
  }
  
  const vector<Photons>& photons(){
    if( Photon_.size() > 0) return Photon_;
    loadPhotons();
  	Photon_.reserve(nPhoton);
    for(size_t idx = 0; idx < nPhoton; ++idx ){
      Photons obj;
      obj.setenergyErr(Photon_energyErr_[idx]);
      obj.sethoe(Photon_hoe_[idx]);
      obj.setmvaID(Photon_mvaID_[idx]);
      obj.setpfRelIso03_all(Photon_pfRelIso03_all_[idx]);
      obj.setpfRelIso03_chg(Photon_pfRelIso03_chg_[idx]);
      obj.setr9(Photon_r9_[idx]);
      obj.setsieie(Photon_sieie_[idx]);
      obj.setcharge(Photon_charge_[idx]);
      obj.setcutBasedBitmap(Photon_cutBasedBitmap_[idx]);
      obj.setelectronIdx(Photon_electronIdx_[idx]);
      obj.setjetIdx(Photon_jetIdx_[idx]);
      obj.setpdgId(Photon_pdgId_[idx]);
      obj.setvidNestedWPBitmap(Photon_vidNestedWPBitmap_[idx]);
      obj.setelectronVeto(Photon_electronVeto_[idx]);
      obj.setisScEtaEB(Photon_isScEtaEB_[idx]);
      obj.setisScEtaEE(Photon_isScEtaEE_[idx]);
      obj.setmvaID_WP80(Photon_mvaID_WP80_[idx]);
      obj.setmvaID_WP90(Photon_mvaID_WP90_[idx]);
      obj.setpixelSeed(Photon_pixelSeed_[idx]);
      obj.setgenPartIdx(Photon_genPartIdx_[idx]);
      obj.setgenPartFlav(Photon_genPartFlav_[idx]);
      obj.setcleanmask(Photon_cleanmask_[idx]);
      obj.setLorentzVector(Photon_pt_[idx], Photon_eta_[idx], Photon_phi_[idx], Photon_mass_[idx]);
      Photon_.push_back( obj );
    }
    return Photon_;
  }
  
  const vector<Genjets>& genjets(){
    if( GenJet_.size() > 0) return GenJet_;
    loadGenjets();
  	GenJet_.reserve(nGenJet);
    for(size_t idx = 0; idx < nGenJet; ++idx ){
      Genjets obj;
      obj.setpartonFlavour(GenJetAK8_partonFlavour_[idx]);
      obj.sethadronFlavour(GenJetAK8_hadronFlavour_[idx]);
      obj.setpartonFlavour(GenJet_partonFlavour_[idx]);
      obj.sethadronFlavour(GenJet_hadronFlavour_[idx]);
      obj.setLorentzVector(GenJet_pt_[idx], GenJet_eta_[idx], GenJet_phi_[idx], GenJet_mass_[idx]);
      GenJet_.push_back( obj );
    }
    return GenJet_;
  }
  
  const Rawmet rawmet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadRawmet();
  
    Rawmet obj;
    obj.setsumEt(RawMET_sumEt_);
    
  
    return obj;
  }
  
  const Btagweight btagweight(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadBtagweight();
  
    Btagweight obj;
    obj.setCSVV2(btagWeight_CSVV2_);
    obj.setDeepCSVB(btagWeight_DeepCSVB_);
    
  
    return obj;
  }
  
  const Softactivityjet softactivityjet(){
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
  
  const L1Simulation l1simulation(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadL1Simulation();
  
    L1Simulation obj;
    obj.setstep(L1simulation_step_);
    
  
    return obj;
  }
  
  const vector<Genparts>& genparts(){
    if( GenPart_.size() > 0) return GenPart_;
    loadGenparts();
  	GenPart_.reserve(nGenPart);
    for(size_t idx = 0; idx < nGenPart; ++idx ){
      Genparts obj;
      obj.setgenPartIdxMother(GenPart_genPartIdxMother_[idx]);
      obj.setpdgId(GenPart_pdgId_[idx]);
      obj.setstatus(GenPart_status_[idx]);
      obj.setstatusFlags(GenPart_statusFlags_[idx]);
      obj.setLorentzVector(GenPart_pt_[idx], GenPart_eta_[idx], GenPart_phi_[idx], GenPart_mass_[idx]);
      GenPart_.push_back( obj );
    }
    return GenPart_;
  }
  
  const Lhe lhe(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadLhe();
  
    Lhe obj;
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
    for( size_t idx = 0; idx < nLHEPdfWeight; ++idx ){
      LHEPdfWeight.push_back(LHEPdfWeight[idx]);
    }
    obj.setLHEPdfWeight(LHEPdfWeight);
    for( size_t idx = 0; idx < nLHEScaleWeight; ++idx ){
      LHEScaleWeight.push_back(LHEScaleWeight[idx]);
    }
    obj.setLHEScaleWeight(LHEScaleWeight);
  
    return obj;
  }
  
  const Tkmet tkmet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadTkmet();
  
    Tkmet obj;
    obj.setsumEt(TkMET_sumEt_);
    
  
    return obj;
  }
  
  const vector<Taus>& taus(){
    if( Tau_.size() > 0) return Tau_;
    loadTaus();
  	Tau_.reserve(nTau);
    for(size_t idx = 0; idx < nTau; ++idx ){
      Taus obj;
      obj.setchargedIso(Tau_chargedIso_[idx]);
      obj.setdxy(Tau_dxy_[idx]);
      obj.setdz(Tau_dz_[idx]);
      obj.setleadTkDeltaEta(Tau_leadTkDeltaEta_[idx]);
      obj.setleadTkDeltaPhi(Tau_leadTkDeltaPhi_[idx]);
      obj.setleadTkPtOverTauPt(Tau_leadTkPtOverTauPt_[idx]);
      obj.setneutralIso(Tau_neutralIso_[idx]);
      obj.setphotonsOutsideSignalCone(Tau_photonsOutsideSignalCone_[idx]);
      obj.setpuCorr(Tau_puCorr_[idx]);
      obj.setrawAntiEle(Tau_rawAntiEle_[idx]);
      obj.setrawIso(Tau_rawIso_[idx]);
      obj.setrawIsodR03(Tau_rawIsodR03_[idx]);
      obj.setcharge(Tau_charge_[idx]);
      obj.setdecayMode(Tau_decayMode_[idx]);
      obj.setjetIdx(Tau_jetIdx_[idx]);
      obj.setrawAntiEleCat(Tau_rawAntiEleCat_[idx]);
      obj.setidAntiEle(Tau_idAntiEle_[idx]);
      obj.setidAntiMu(Tau_idAntiMu_[idx]);
      obj.setidDecayMode(Tau_idDecayMode_[idx]);
      obj.setidDecayModeNewDMs(Tau_idDecayModeNewDMs_[idx]);
      obj.setcleanmask(Tau_cleanmask_[idx]);
      obj.setgenPartIdx(Tau_genPartIdx_[idx]);
      obj.setgenPartFlav(Tau_genPartFlav_[idx]);
      obj.setLorentzVector(Tau_pt_[idx], Tau_eta_[idx], Tau_phi_[idx], Tau_mass_[idx]);
      Tau_.push_back( obj );
    }
    return Tau_;
  }
  
  const Flag flag(){
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
    obj.setecalBadCalibFilter(Flag_ecalBadCalibFilter_);
    obj.setgoodVertices(Flag_goodVertices_);
    obj.seteeBadScFilter(Flag_eeBadScFilter_);
    obj.setecalLaserCorrFilter(Flag_ecalLaserCorrFilter_);
    obj.settrkPOGFilters(Flag_trkPOGFilters_);
    obj.setchargedHadronTrackResolutionFilter(Flag_chargedHadronTrackResolutionFilter_);
    obj.setmuonBadTrackFilter(Flag_muonBadTrackFilter_);
    obj.setBadChargedCandidateFilter(Flag_BadChargedCandidateFilter_);
    obj.setBadPFMuonFilter(Flag_BadPFMuonFilter_);
    obj.setBadChargedCandidateSummer16Filter(Flag_BadChargedCandidateSummer16Filter_);
    obj.setBadPFMuonSummer16Filter(Flag_BadPFMuonSummer16Filter_);
    obj.settrkPOG_manystripclus53X(Flag_trkPOG_manystripclus53X_);
    obj.settrkPOG_toomanystripclus53X(Flag_trkPOG_toomanystripclus53X_);
    obj.settrkPOG_logErrorTooManyClusters(Flag_trkPOG_logErrorTooManyClusters_);
    obj.setMETFilters(Flag_METFilters_);
    
  
    return obj;
  }
  
  const Puppimet puppimet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadPuppimet();
  
    Puppimet obj;
    obj.setsumEt(PuppiMET_sumEt_);
    
  
    return obj;
  }
  
  const vector<Muons>& muons(){
    if( Muon_.size() > 0) return Muon_;
    loadMuons();
  	Muon_.reserve(nMuon);
    for(size_t idx = 0; idx < nMuon; ++idx ){
      Muons obj;
      obj.setdxy(Muon_dxy_[idx]);
      obj.setdxyErr(Muon_dxyErr_[idx]);
      obj.setdz(Muon_dz_[idx]);
      obj.setdzErr(Muon_dzErr_[idx]);
      obj.setip3d(Muon_ip3d_[idx]);
      obj.setminiPFRelIso_all(Muon_miniPFRelIso_all_[idx]);
      obj.setminiPFRelIso_chg(Muon_miniPFRelIso_chg_[idx]);
      obj.setpfRelIso03_all(Muon_pfRelIso03_all_[idx]);
      obj.setpfRelIso03_chg(Muon_pfRelIso03_chg_[idx]);
      obj.setpfRelIso04_all(Muon_pfRelIso04_all_[idx]);
      obj.setptErr(Muon_ptErr_[idx]);
      obj.setsegmentComp(Muon_segmentComp_[idx]);
      obj.setsip3d(Muon_sip3d_[idx]);
      obj.setmvaTTH(Muon_mvaTTH_[idx]);
      obj.setcharge(Muon_charge_[idx]);
      obj.setjetIdx(Muon_jetIdx_[idx]);
      obj.setnStations(Muon_nStations_[idx]);
      obj.setnTrackerLayers(Muon_nTrackerLayers_[idx]);
      obj.setpdgId(Muon_pdgId_[idx]);
      obj.settightCharge(Muon_tightCharge_[idx]);
      obj.sethighPtId(Muon_highPtId_[idx]);
      obj.setisPFcand(Muon_isPFcand_[idx]);
      obj.setmediumId(Muon_mediumId_[idx]);
      obj.setsoftId(Muon_softId_[idx]);
      obj.settightId(Muon_tightId_[idx]);
      obj.setgenPartIdx(Muon_genPartIdx_[idx]);
      obj.setgenPartFlav(Muon_genPartFlav_[idx]);
      obj.setcleanmask(Muon_cleanmask_[idx]);
      obj.setLorentzVector(Muon_pt_[idx], Muon_eta_[idx], Muon_phi_[idx], Muon_mass_[idx]);
      Muon_.push_back( obj );
    }
    return Muon_;
  }
  
  const L1Reco l1reco(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadL1Reco();
  
    L1Reco obj;
    obj.setstep(L1Reco_step_);
    
  
    return obj;
  }
  
  const vector<Otherpvs>& otherpvs(){
    if( OtherPV_.size() > 0) return OtherPV_;
    loadOtherpvs();
  	OtherPV_.reserve(nOtherPV);
    for(size_t idx = 0; idx < nOtherPV; ++idx ){
      Otherpvs obj;
      obj.setz(OtherPV_z_[idx]);
      
      OtherPV_.push_back( obj );
    }
    return OtherPV_;
  }
  
  const Hlt hlt(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadHlt();
  
    Hlt obj;
    obj.setHLTriggerFirstPath(HLTriggerFirstPath);
    obj.setAK8PFJet360_TrimMass30(HLT_AK8PFJet360_TrimMass30_);
    obj.setAK8PFJet380_TrimMass30(HLT_AK8PFJet380_TrimMass30_);
    obj.setAK8PFJet400_TrimMass30(HLT_AK8PFJet400_TrimMass30_);
    obj.setAK8PFJet420_TrimMass30(HLT_AK8PFJet420_TrimMass30_);
    obj.setAK8PFHT750_TrimMass50(HLT_AK8PFHT750_TrimMass50_);
    obj.setAK8PFHT800_TrimMass50(HLT_AK8PFHT800_TrimMass50_);
    obj.setAK8PFHT850_TrimMass50(HLT_AK8PFHT850_TrimMass50_);
    obj.setAK8PFHT900_TrimMass50(HLT_AK8PFHT900_TrimMass50_);
    obj.setCaloJet500_NoJetID(HLT_CaloJet500_NoJetID_);
    obj.setCaloJet550_NoJetID(HLT_CaloJet550_NoJetID_);
    obj.setDoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL(HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_);
    obj.setDoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon(HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_);
    obj.setTrimuon5_3p5_2_Upsilon_Muon(HLT_Trimuon5_3p5_2_Upsilon_Muon_);
    obj.setTrimuonOpen_5_3p5_2_Upsilon_Muon(HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_);
    obj.setDoubleEle25_CaloIdL_MW(HLT_DoubleEle25_CaloIdL_MW_);
    obj.setDoubleEle27_CaloIdL_MW(HLT_DoubleEle27_CaloIdL_MW_);
    obj.setDoubleEle33_CaloIdL_MW(HLT_DoubleEle33_CaloIdL_MW_);
    obj.setDoubleEle24_eta2p1_WPTight_Gsf(HLT_DoubleEle24_eta2p1_WPTight_Gsf_);
    obj.setDoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350(HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_);
    obj.setDoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350(HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_);
    obj.setEle27_Ele37_CaloIdL_MW(HLT_Ele27_Ele37_CaloIdL_MW_);
    obj.setMu27_Ele37_CaloIdL_MW(HLT_Mu27_Ele37_CaloIdL_MW_);
    obj.setMu37_Ele27_CaloIdL_MW(HLT_Mu37_Ele27_CaloIdL_MW_);
    obj.setMu37_TkMu27(HLT_Mu37_TkMu27_);
    obj.setDoubleMu4_3_Bs(HLT_DoubleMu4_3_Bs_);
    obj.setDoubleMu4_3_Jpsi(HLT_DoubleMu4_3_Jpsi_);
    obj.setDoubleMu4_JpsiTrk_Displaced(HLT_DoubleMu4_JpsiTrk_Displaced_);
    obj.setDoubleMu4_LowMassNonResonantTrk_Displaced(HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_);
    obj.setDoubleMu3_Trk_Tau3mu(HLT_DoubleMu3_Trk_Tau3mu_);
    obj.setDoubleMu3_TkMu_DsTau3Mu(HLT_DoubleMu3_TkMu_DsTau3Mu_);
    obj.setDoubleMu4_PsiPrimeTrk_Displaced(HLT_DoubleMu4_PsiPrimeTrk_Displaced_);
    obj.setDoubleMu4_Mass3p8_DZ_PFHT350(HLT_DoubleMu4_Mass3p8_DZ_PFHT350_);
    obj.setMu3_PFJet40(HLT_Mu3_PFJet40_);
    obj.setMu7p5_L2Mu2_Jpsi(HLT_Mu7p5_L2Mu2_Jpsi_);
    obj.setMu7p5_L2Mu2_Upsilon(HLT_Mu7p5_L2Mu2_Upsilon_);
    obj.setMu7p5_Track2_Jpsi(HLT_Mu7p5_Track2_Jpsi_);
    obj.setMu7p5_Track3p5_Jpsi(HLT_Mu7p5_Track3p5_Jpsi_);
    obj.setMu7p5_Track7_Jpsi(HLT_Mu7p5_Track7_Jpsi_);
    obj.setMu7p5_Track2_Upsilon(HLT_Mu7p5_Track2_Upsilon_);
    obj.setMu7p5_Track3p5_Upsilon(HLT_Mu7p5_Track3p5_Upsilon_);
    obj.setMu7p5_Track7_Upsilon(HLT_Mu7p5_Track7_Upsilon_);
    obj.setMu3_L1SingleMu5orSingleMu7(HLT_Mu3_L1SingleMu5orSingleMu7_);
    obj.setDoublePhoton33_CaloIdL(HLT_DoublePhoton33_CaloIdL_);
    obj.setDoublePhoton70(HLT_DoublePhoton70_);
    obj.setDoublePhoton85(HLT_DoublePhoton85_);
    obj.setEle20_WPTight_Gsf(HLT_Ele20_WPTight_Gsf_);
    obj.setEle15_WPLoose_Gsf(HLT_Ele15_WPLoose_Gsf_);
    obj.setEle17_WPLoose_Gsf(HLT_Ele17_WPLoose_Gsf_);
    obj.setEle20_WPLoose_Gsf(HLT_Ele20_WPLoose_Gsf_);
    obj.setEle20_eta2p1_WPLoose_Gsf(HLT_Ele20_eta2p1_WPLoose_Gsf_);
    obj.setDiEle27_WPTightCaloOnly_L1DoubleEG(HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_);
    obj.setEle27_WPTight_Gsf(HLT_Ele27_WPTight_Gsf_);
    obj.setEle28_WPTight_Gsf(HLT_Ele28_WPTight_Gsf_);
    obj.setEle30_WPTight_Gsf(HLT_Ele30_WPTight_Gsf_);
    obj.setEle32_WPTight_Gsf(HLT_Ele32_WPTight_Gsf_);
    obj.setEle35_WPTight_Gsf(HLT_Ele35_WPTight_Gsf_);
    obj.setEle35_WPTight_Gsf_L1EGMT(HLT_Ele35_WPTight_Gsf_L1EGMT_);
    obj.setEle38_WPTight_Gsf(HLT_Ele38_WPTight_Gsf_);
    obj.setEle40_WPTight_Gsf(HLT_Ele40_WPTight_Gsf_);
    obj.setEle32_WPTight_Gsf_L1DoubleEG(HLT_Ele32_WPTight_Gsf_L1DoubleEG_);
    obj.setEle24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1(HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_);
    obj.setEle24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1(HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_);
    obj.setEle24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1(HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_);
    obj.setEle24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1(HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_);
    obj.setEle24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1(HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_);
    obj.setEle24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1(HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_);
    obj.setHT450_Beamspot(HLT_HT450_Beamspot_);
    obj.setHT300_Beamspot(HLT_HT300_Beamspot_);
    obj.setZeroBias_Beamspot(HLT_ZeroBias_Beamspot_);
    obj.setIsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1(HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_);
    obj.setIsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1(HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_);
    obj.setIsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1(HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_);
    obj.setIsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1(HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_);
    obj.setIsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1(HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_);
    obj.setIsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1(HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_);
    obj.setIsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1(HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_);
    obj.setIsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1(HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_);
    obj.setIsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1(HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_);
    obj.setIsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1(HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_);
    obj.setIsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1(HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_);
    obj.setIsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1(HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_);
    obj.setIsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1(HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_);
    obj.setIsoMu20(HLT_IsoMu20_);
    obj.setIsoMu24(HLT_IsoMu24_);
    obj.setIsoMu24_eta2p1(HLT_IsoMu24_eta2p1_);
    obj.setIsoMu27(HLT_IsoMu27_);
    obj.setIsoMu30(HLT_IsoMu30_);
    obj.setUncorrectedJetE30_NoBPTX(HLT_UncorrectedJetE30_NoBPTX_);
    obj.setUncorrectedJetE30_NoBPTX3BX(HLT_UncorrectedJetE30_NoBPTX3BX_);
    obj.setUncorrectedJetE60_NoBPTX3BX(HLT_UncorrectedJetE60_NoBPTX3BX_);
    obj.setUncorrectedJetE70_NoBPTX3BX(HLT_UncorrectedJetE70_NoBPTX3BX_);
    obj.setL1SingleMu18(HLT_L1SingleMu18_);
    obj.setL1SingleMu25(HLT_L1SingleMu25_);
    obj.setL2Mu10(HLT_L2Mu10_);
    obj.setL2Mu10_NoVertex_NoBPTX3BX(HLT_L2Mu10_NoVertex_NoBPTX3BX_);
    obj.setL2Mu10_NoVertex_NoBPTX(HLT_L2Mu10_NoVertex_NoBPTX_);
    obj.setL2Mu45_NoVertex_3Sta_NoBPTX3BX(HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_);
    obj.setL2Mu40_NoVertex_3Sta_NoBPTX3BX(HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_);
    obj.setL2Mu50(HLT_L2Mu50_);
    obj.setL2Mu23NoVtx_2Cha(HLT_L2Mu23NoVtx_2Cha_);
    obj.setL2Mu23NoVtx_2Cha_CosmicSeed(HLT_L2Mu23NoVtx_2Cha_CosmicSeed_);
    obj.setDoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4(HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_);
    obj.setDoubleL2Mu30NoVtx_2Cha_Eta2p4(HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_);
    obj.setDoubleL2Mu50(HLT_DoubleL2Mu50_);
    obj.setDoubleL2Mu23NoVtx_2Cha_CosmicSeed(HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_);
    obj.setDoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched(HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_);
    obj.setDoubleL2Mu25NoVtx_2Cha_CosmicSeed(HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_);
    obj.setDoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched(HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_);
    obj.setDoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4(HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_);
    obj.setDoubleL2Mu23NoVtx_2Cha(HLT_DoubleL2Mu23NoVtx_2Cha_);
    obj.setDoubleL2Mu23NoVtx_2Cha_NoL2Matched(HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_);
    obj.setDoubleL2Mu25NoVtx_2Cha(HLT_DoubleL2Mu25NoVtx_2Cha_);
    obj.setDoubleL2Mu25NoVtx_2Cha_NoL2Matched(HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_);
    obj.setDoubleL2Mu25NoVtx_2Cha_Eta2p4(HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_);
    obj.setMu17_TrkIsoVVL_Mu8_TrkIsoVVL(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_);
    obj.setMu19_TrkIsoVVL_Mu9_TrkIsoVVL(HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_);
    obj.setMu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_);
    obj.setMu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ(HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_);
    obj.setMu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_);
    obj.setMu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8(HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_);
    obj.setMu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_);
    obj.setMu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8(HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_);
    obj.setMu25_TkMu0_Onia(HLT_Mu25_TkMu0_Onia_);
    obj.setMu30_TkMu0_Psi(HLT_Mu30_TkMu0_Psi_);
    obj.setMu30_TkMu0_Upsilon(HLT_Mu30_TkMu0_Upsilon_);
    obj.setMu20_TkMu0_Phi(HLT_Mu20_TkMu0_Phi_);
    obj.setMu25_TkMu0_Phi(HLT_Mu25_TkMu0_Phi_);
    obj.setMu12(HLT_Mu12_);
    obj.setMu15(HLT_Mu15_);
    obj.setMu20(HLT_Mu20_);
    obj.setMu27(HLT_Mu27_);
    obj.setMu50(HLT_Mu50_);
    obj.setMu55(HLT_Mu55_);
    obj.setOldMu100(HLT_OldMu100_);
    obj.setTkMu100(HLT_TkMu100_);
    obj.setDiPFJetAve40(HLT_DiPFJetAve40_);
    obj.setDiPFJetAve60(HLT_DiPFJetAve60_);
    obj.setDiPFJetAve80(HLT_DiPFJetAve80_);
    obj.setDiPFJetAve140(HLT_DiPFJetAve140_);
    obj.setDiPFJetAve200(HLT_DiPFJetAve200_);
    obj.setDiPFJetAve260(HLT_DiPFJetAve260_);
    obj.setDiPFJetAve320(HLT_DiPFJetAve320_);
    obj.setDiPFJetAve400(HLT_DiPFJetAve400_);
    obj.setDiPFJetAve500(HLT_DiPFJetAve500_);
    obj.setDiPFJetAve60_HFJEC(HLT_DiPFJetAve60_HFJEC_);
    obj.setDiPFJetAve80_HFJEC(HLT_DiPFJetAve80_HFJEC_);
    obj.setDiPFJetAve100_HFJEC(HLT_DiPFJetAve100_HFJEC_);
    obj.setDiPFJetAve160_HFJEC(HLT_DiPFJetAve160_HFJEC_);
    obj.setDiPFJetAve220_HFJEC(HLT_DiPFJetAve220_HFJEC_);
    obj.setDiPFJetAve300_HFJEC(HLT_DiPFJetAve300_HFJEC_);
    obj.setAK8PFJet15(HLT_AK8PFJet15_);
    obj.setAK8PFJet25(HLT_AK8PFJet25_);
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
    obj.setAK8PFJet550(HLT_AK8PFJet550_);
    obj.setPFJet15(HLT_PFJet15_);
    obj.setPFJet25(HLT_PFJet25_);
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
    obj.setPFJet550(HLT_PFJet550_);
    obj.setPFJetFwd15(HLT_PFJetFwd15_);
    obj.setPFJetFwd25(HLT_PFJetFwd25_);
    obj.setPFJetFwd40(HLT_PFJetFwd40_);
    obj.setPFJetFwd60(HLT_PFJetFwd60_);
    obj.setPFJetFwd80(HLT_PFJetFwd80_);
    obj.setPFJetFwd140(HLT_PFJetFwd140_);
    obj.setPFJetFwd200(HLT_PFJetFwd200_);
    obj.setPFJetFwd260(HLT_PFJetFwd260_);
    obj.setPFJetFwd320(HLT_PFJetFwd320_);
    obj.setPFJetFwd400(HLT_PFJetFwd400_);
    obj.setPFJetFwd450(HLT_PFJetFwd450_);
    obj.setPFJetFwd500(HLT_PFJetFwd500_);
    obj.setAK8PFJetFwd15(HLT_AK8PFJetFwd15_);
    obj.setAK8PFJetFwd25(HLT_AK8PFJetFwd25_);
    obj.setAK8PFJetFwd40(HLT_AK8PFJetFwd40_);
    obj.setAK8PFJetFwd60(HLT_AK8PFJetFwd60_);
    obj.setAK8PFJetFwd80(HLT_AK8PFJetFwd80_);
    obj.setAK8PFJetFwd140(HLT_AK8PFJetFwd140_);
    obj.setAK8PFJetFwd200(HLT_AK8PFJetFwd200_);
    obj.setAK8PFJetFwd260(HLT_AK8PFJetFwd260_);
    obj.setAK8PFJetFwd320(HLT_AK8PFJetFwd320_);
    obj.setAK8PFJetFwd400(HLT_AK8PFJetFwd400_);
    obj.setAK8PFJetFwd450(HLT_AK8PFJetFwd450_);
    obj.setAK8PFJetFwd500(HLT_AK8PFJetFwd500_);
    obj.setPFHT180(HLT_PFHT180_);
    obj.setPFHT250(HLT_PFHT250_);
    obj.setPFHT370(HLT_PFHT370_);
    obj.setPFHT430(HLT_PFHT430_);
    obj.setPFHT510(HLT_PFHT510_);
    obj.setPFHT590(HLT_PFHT590_);
    obj.setPFHT680(HLT_PFHT680_);
    obj.setPFHT780(HLT_PFHT780_);
    obj.setPFHT890(HLT_PFHT890_);
    obj.setPFHT1050(HLT_PFHT1050_);
    obj.setPFHT500_PFMET100_PFMHT100_IDTight(HLT_PFHT500_PFMET100_PFMHT100_IDTight_);
    obj.setPFHT500_PFMET110_PFMHT110_IDTight(HLT_PFHT500_PFMET110_PFMHT110_IDTight_);
    obj.setPFHT700_PFMET85_PFMHT85_IDTight(HLT_PFHT700_PFMET85_PFMHT85_IDTight_);
    obj.setPFHT700_PFMET95_PFMHT95_IDTight(HLT_PFHT700_PFMET95_PFMHT95_IDTight_);
    obj.setPFHT800_PFMET75_PFMHT75_IDTight(HLT_PFHT800_PFMET75_PFMHT75_IDTight_);
    obj.setPFHT800_PFMET85_PFMHT85_IDTight(HLT_PFHT800_PFMET85_PFMHT85_IDTight_);
    obj.setPFMET110_PFMHT110_IDTight(HLT_PFMET110_PFMHT110_IDTight_);
    obj.setPFMET120_PFMHT120_IDTight(HLT_PFMET120_PFMHT120_IDTight_);
    obj.setPFMET130_PFMHT130_IDTight(HLT_PFMET130_PFMHT130_IDTight_);
    obj.setPFMET140_PFMHT140_IDTight(HLT_PFMET140_PFMHT140_IDTight_);
    obj.setPFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1(HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_);
    obj.setPFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1(HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_);
    obj.setPFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1(HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_);
    obj.setPFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1(HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_);
    obj.setPFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1(HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_);
    obj.setPFMET120_PFMHT120_IDTight_PFHT60(HLT_PFMET120_PFMHT120_IDTight_PFHT60_);
    obj.setPFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_);
    obj.setPFMETTypeOne120_PFMHT120_IDTight_PFHT60(HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_);
    obj.setPFMETTypeOne110_PFMHT110_IDTight(HLT_PFMETTypeOne110_PFMHT110_IDTight_);
    obj.setPFMETTypeOne120_PFMHT120_IDTight(HLT_PFMETTypeOne120_PFMHT120_IDTight_);
    obj.setPFMETTypeOne130_PFMHT130_IDTight(HLT_PFMETTypeOne130_PFMHT130_IDTight_);
    obj.setPFMETTypeOne140_PFMHT140_IDTight(HLT_PFMETTypeOne140_PFMHT140_IDTight_);
    obj.setPFMETNoMu110_PFMHTNoMu110_IDTight(HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_);
    obj.setPFMETNoMu120_PFMHTNoMu120_IDTight(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_);
    obj.setPFMETNoMu130_PFMHTNoMu130_IDTight(HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_);
    obj.setPFMETNoMu140_PFMHTNoMu140_IDTight(HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_);
    obj.setMonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight(HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_);
    obj.setMonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight(HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_);
    obj.setMonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight(HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_);
    obj.setMonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight(HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_);
    obj.setL1ETMHadSeeds(HLT_L1ETMHadSeeds_);
    obj.setCaloMHT90(HLT_CaloMHT90_);
    obj.setCaloMET80_NotCleaned(HLT_CaloMET80_NotCleaned_);
    obj.setCaloMET90_NotCleaned(HLT_CaloMET90_NotCleaned_);
    obj.setCaloMET100_NotCleaned(HLT_CaloMET100_NotCleaned_);
    obj.setCaloMET110_NotCleaned(HLT_CaloMET110_NotCleaned_);
    obj.setCaloMET250_NotCleaned(HLT_CaloMET250_NotCleaned_);
    obj.setCaloMET70_HBHECleaned(HLT_CaloMET70_HBHECleaned_);
    obj.setCaloMET80_HBHECleaned(HLT_CaloMET80_HBHECleaned_);
    obj.setCaloMET90_HBHECleaned(HLT_CaloMET90_HBHECleaned_);
    obj.setCaloMET100_HBHECleaned(HLT_CaloMET100_HBHECleaned_);
    obj.setCaloMET250_HBHECleaned(HLT_CaloMET250_HBHECleaned_);
    obj.setCaloMET300_HBHECleaned(HLT_CaloMET300_HBHECleaned_);
    obj.setCaloMET350_HBHECleaned(HLT_CaloMET350_HBHECleaned_);
    obj.setPFMET200_NotCleaned(HLT_PFMET200_NotCleaned_);
    obj.setPFMET200_HBHECleaned(HLT_PFMET200_HBHECleaned_);
    obj.setPFMET250_HBHECleaned(HLT_PFMET250_HBHECleaned_);
    obj.setPFMET300_HBHECleaned(HLT_PFMET300_HBHECleaned_);
    obj.setPFMET200_HBHE_BeamHaloCleaned(HLT_PFMET200_HBHE_BeamHaloCleaned_);
    obj.setPFMETTypeOne200_HBHE_BeamHaloCleaned(HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_);
    obj.setMET105_IsoTrk50(HLT_MET105_IsoTrk50_);
    obj.setMET120_IsoTrk50(HLT_MET120_IsoTrk50_);
    obj.setSingleJet30_Mu12_SinglePFJet40(HLT_SingleJet30_Mu12_SinglePFJet40_);
    obj.setMu12_DoublePFJets40_CaloBTagDeepCSV_p71(HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_);
    obj.setMu12_DoublePFJets100_CaloBTagDeepCSV_p71(HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_);
    obj.setMu12_DoublePFJets200_CaloBTagDeepCSV_p71(HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_);
    obj.setMu12_DoublePFJets350_CaloBTagDeepCSV_p71(HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_);
    obj.setMu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
    obj.setMu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
    obj.setMu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
    obj.setDoublePFJets40_CaloBTagDeepCSV_p71(HLT_DoublePFJets40_CaloBTagDeepCSV_p71_);
    obj.setDoublePFJets100_CaloBTagDeepCSV_p71(HLT_DoublePFJets100_CaloBTagDeepCSV_p71_);
    obj.setDoublePFJets200_CaloBTagDeepCSV_p71(HLT_DoublePFJets200_CaloBTagDeepCSV_p71_);
    obj.setDoublePFJets350_CaloBTagDeepCSV_p71(HLT_DoublePFJets350_CaloBTagDeepCSV_p71_);
    obj.setDoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
    obj.setDoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71(HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_);
    obj.setPhoton300_NoHE(HLT_Photon300_NoHE_);
    obj.setMu8_TrkIsoVVL(HLT_Mu8_TrkIsoVVL_);
    obj.setMu8_DiEle12_CaloIdL_TrackIdL_DZ(HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_);
    obj.setMu8_DiEle12_CaloIdL_TrackIdL(HLT_Mu8_DiEle12_CaloIdL_TrackIdL_);
    obj.setMu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ(HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_);
    obj.setMu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350(HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_);
    obj.setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_);
    obj.setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_);
    obj.setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_);
    obj.setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_);
    obj.setMu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu17_TrkIsoVVL(HLT_Mu17_TrkIsoVVL_);
    obj.setMu19_TrkIsoVVL(HLT_Mu19_TrkIsoVVL_);
    obj.setBTagMu_AK4DiJet20_Mu5(HLT_BTagMu_AK4DiJet20_Mu5_);
    obj.setBTagMu_AK4DiJet40_Mu5(HLT_BTagMu_AK4DiJet40_Mu5_);
    obj.setBTagMu_AK4DiJet70_Mu5(HLT_BTagMu_AK4DiJet70_Mu5_);
    obj.setBTagMu_AK4DiJet110_Mu5(HLT_BTagMu_AK4DiJet110_Mu5_);
    obj.setBTagMu_AK4DiJet170_Mu5(HLT_BTagMu_AK4DiJet170_Mu5_);
    obj.setBTagMu_AK4Jet300_Mu5(HLT_BTagMu_AK4Jet300_Mu5_);
    obj.setBTagMu_AK8DiJet170_Mu5(HLT_BTagMu_AK8DiJet170_Mu5_);
    obj.setBTagMu_AK8Jet170_DoubleMu5(HLT_BTagMu_AK8Jet170_DoubleMu5_);
    obj.setBTagMu_AK8Jet300_Mu5(HLT_BTagMu_AK8Jet300_Mu5_);
    obj.setBTagMu_AK4DiJet20_Mu5_noalgo(HLT_BTagMu_AK4DiJet20_Mu5_noalgo_);
    obj.setBTagMu_AK4DiJet40_Mu5_noalgo(HLT_BTagMu_AK4DiJet40_Mu5_noalgo_);
    obj.setBTagMu_AK4DiJet70_Mu5_noalgo(HLT_BTagMu_AK4DiJet70_Mu5_noalgo_);
    obj.setBTagMu_AK4DiJet110_Mu5_noalgo(HLT_BTagMu_AK4DiJet110_Mu5_noalgo_);
    obj.setBTagMu_AK4DiJet170_Mu5_noalgo(HLT_BTagMu_AK4DiJet170_Mu5_noalgo_);
    obj.setBTagMu_AK4Jet300_Mu5_noalgo(HLT_BTagMu_AK4Jet300_Mu5_noalgo_);
    obj.setBTagMu_AK8DiJet170_Mu5_noalgo(HLT_BTagMu_AK8DiJet170_Mu5_noalgo_);
    obj.setBTagMu_AK8Jet170_DoubleMu5_noalgo(HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_);
    obj.setBTagMu_AK8Jet300_Mu5_noalgo(HLT_BTagMu_AK8Jet300_Mu5_noalgo_);
    obj.setEle15_Ele8_CaloIdL_TrackIdL_IsoVL(HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_);
    obj.setEle23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setEle23_Ele12_CaloIdL_TrackIdL_IsoVL(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setMu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL(HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_);
    obj.setMu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ(HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_);
    obj.setMu12_DoublePhoton20(HLT_Mu12_DoublePhoton20_);
    obj.setTriplePhoton_20_20_20_CaloIdLV2(HLT_TriplePhoton_20_20_20_CaloIdLV2_);
    obj.setTriplePhoton_20_20_20_CaloIdLV2_R9IdVL(HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_);
    obj.setTriplePhoton_30_30_10_CaloIdLV2(HLT_TriplePhoton_30_30_10_CaloIdLV2_);
    obj.setTriplePhoton_30_30_10_CaloIdLV2_R9IdVL(HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_);
    obj.setTriplePhoton_35_35_5_CaloIdLV2_R9IdVL(HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_);
    obj.setPhoton20(HLT_Photon20_);
    obj.setPhoton33(HLT_Photon33_);
    obj.setPhoton50(HLT_Photon50_);
    obj.setPhoton75(HLT_Photon75_);
    obj.setPhoton90(HLT_Photon90_);
    obj.setPhoton120(HLT_Photon120_);
    obj.setPhoton150(HLT_Photon150_);
    obj.setPhoton175(HLT_Photon175_);
    obj.setPhoton200(HLT_Photon200_);
    obj.setPhoton100EB_TightID_TightIso(HLT_Photon100EB_TightID_TightIso_);
    obj.setPhoton110EB_TightID_TightIso(HLT_Photon110EB_TightID_TightIso_);
    obj.setPhoton120EB_TightID_TightIso(HLT_Photon120EB_TightID_TightIso_);
    obj.setPhoton100EBHE10(HLT_Photon100EBHE10_);
    obj.setPhoton100EEHE10(HLT_Photon100EEHE10_);
    obj.setPhoton100EE_TightID_TightIso(HLT_Photon100EE_TightID_TightIso_);
    obj.setPhoton50_R9Id90_HE10_IsoM(HLT_Photon50_R9Id90_HE10_IsoM_);
    obj.setPhoton75_R9Id90_HE10_IsoM(HLT_Photon75_R9Id90_HE10_IsoM_);
    obj.setPhoton75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3(HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_);
    obj.setPhoton75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3(HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_);
    obj.setPhoton90_R9Id90_HE10_IsoM(HLT_Photon90_R9Id90_HE10_IsoM_);
    obj.setPhoton120_R9Id90_HE10_IsoM(HLT_Photon120_R9Id90_HE10_IsoM_);
    obj.setPhoton165_R9Id90_HE10_IsoM(HLT_Photon165_R9Id90_HE10_IsoM_);
    obj.setPhoton90_CaloIdL_PFHT700(HLT_Photon90_CaloIdL_PFHT700_);
    obj.setDiphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90(HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_);
    obj.setDiphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95(HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_);
    obj.setDiphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55(HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_);
    obj.setDiphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55(HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_);
    obj.setPhoton35_TwoProngs35(HLT_Photon35_TwoProngs35_);
    obj.setIsoMu24_TwoProngs35(HLT_IsoMu24_TwoProngs35_);
    obj.setDimuon0_Jpsi_L1_NoOS(HLT_Dimuon0_Jpsi_L1_NoOS_);
    obj.setDimuon0_Jpsi_NoVertexing_NoOS(HLT_Dimuon0_Jpsi_NoVertexing_NoOS_);
    obj.setDimuon0_Jpsi(HLT_Dimuon0_Jpsi_);
    obj.setDimuon0_Jpsi_NoVertexing(HLT_Dimuon0_Jpsi_NoVertexing_);
    obj.setDimuon0_Jpsi_L1_4R_0er1p5R(HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_);
    obj.setDimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R(HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_);
    obj.setDimuon0_Jpsi3p5_Muon2(HLT_Dimuon0_Jpsi3p5_Muon2_);
    obj.setDimuon0_Upsilon_L1_4p5(HLT_Dimuon0_Upsilon_L1_4p5_);
    obj.setDimuon0_Upsilon_L1_5(HLT_Dimuon0_Upsilon_L1_5_);
    obj.setDimuon0_Upsilon_L1_4p5NoOS(HLT_Dimuon0_Upsilon_L1_4p5NoOS_);
    obj.setDimuon0_Upsilon_L1_4p5er2p0(HLT_Dimuon0_Upsilon_L1_4p5er2p0_);
    obj.setDimuon0_Upsilon_L1_4p5er2p0M(HLT_Dimuon0_Upsilon_L1_4p5er2p0M_);
    obj.setDimuon0_Upsilon_NoVertexing(HLT_Dimuon0_Upsilon_NoVertexing_);
    obj.setDimuon0_Upsilon_L1_5M(HLT_Dimuon0_Upsilon_L1_5M_);
    obj.setDimuon0_LowMass_L1_0er1p5R(HLT_Dimuon0_LowMass_L1_0er1p5R_);
    obj.setDimuon0_LowMass_L1_0er1p5(HLT_Dimuon0_LowMass_L1_0er1p5_);
    obj.setDimuon0_LowMass(HLT_Dimuon0_LowMass_);
    obj.setDimuon0_LowMass_L1_4(HLT_Dimuon0_LowMass_L1_4_);
    obj.setDimuon0_LowMass_L1_4R(HLT_Dimuon0_LowMass_L1_4R_);
    obj.setDimuon0_LowMass_L1_TM530(HLT_Dimuon0_LowMass_L1_TM530_);
    obj.setDimuon0_Upsilon_Muon_L1_TM0(HLT_Dimuon0_Upsilon_Muon_L1_TM0_);
    obj.setDimuon0_Upsilon_Muon_NoL1Mass(HLT_Dimuon0_Upsilon_Muon_NoL1Mass_);
    obj.setTripleMu_5_3_3_Mass3p8_DZ(HLT_TripleMu_5_3_3_Mass3p8_DZ_);
    obj.setTripleMu_10_5_5_DZ(HLT_TripleMu_10_5_5_DZ_);
    obj.setTripleMu_12_10_5(HLT_TripleMu_12_10_5_);
    obj.setTau3Mu_Mu7_Mu1_TkMu1_Tau15(HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_);
    obj.setTau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1(HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_);
    obj.setTau3Mu_Mu7_Mu1_TkMu1_IsoTau15(HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_);
    obj.setTau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1(HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_);
    obj.setDoubleMu3_DZ_PFMET50_PFMHT60(HLT_DoubleMu3_DZ_PFMET50_PFMHT60_);
    obj.setDoubleMu3_DZ_PFMET70_PFMHT70(HLT_DoubleMu3_DZ_PFMET70_PFMHT70_);
    obj.setDoubleMu3_DZ_PFMET90_PFMHT90(HLT_DoubleMu3_DZ_PFMET90_PFMHT90_);
    obj.setDoubleMu3_Trk_Tau3mu_NoL1Mass(HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_);
    obj.setDoubleMu4_Jpsi_Displaced(HLT_DoubleMu4_Jpsi_Displaced_);
    obj.setDoubleMu4_Jpsi_NoVertexing(HLT_DoubleMu4_Jpsi_NoVertexing_);
    obj.setDoubleMu4_JpsiTrkTrk_Displaced(HLT_DoubleMu4_JpsiTrkTrk_Displaced_);
    obj.setDoubleMu43NoFiltersNoVtx(HLT_DoubleMu43NoFiltersNoVtx_);
    obj.setDoubleMu48NoFiltersNoVtx(HLT_DoubleMu48NoFiltersNoVtx_);
    obj.setMu43NoFiltersNoVtx_Photon43_CaloIdL(HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_);
    obj.setMu48NoFiltersNoVtx_Photon48_CaloIdL(HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_);
    obj.setMu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL(HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_);
    obj.setMu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL(HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_);
    obj.setDoubleMu33NoFiltersNoVtxDisplaced(HLT_DoubleMu33NoFiltersNoVtxDisplaced_);
    obj.setDoubleMu40NoFiltersNoVtxDisplaced(HLT_DoubleMu40NoFiltersNoVtxDisplaced_);
    obj.setDoubleMu20_7_Mass0to30_L1_DM4(HLT_DoubleMu20_7_Mass0to30_L1_DM4_);
    obj.setDoubleMu20_7_Mass0to30_L1_DM4EG(HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_);
    obj.setHT425(HLT_HT425_);
    obj.setHT430_DisplacedDijet40_DisplacedTrack(HLT_HT430_DisplacedDijet40_DisplacedTrack_);
    obj.setHT500_DisplacedDijet40_DisplacedTrack(HLT_HT500_DisplacedDijet40_DisplacedTrack_);
    obj.setHT430_DisplacedDijet60_DisplacedTrack(HLT_HT430_DisplacedDijet60_DisplacedTrack_);
    obj.setHT400_DisplacedDijet40_DisplacedTrack(HLT_HT400_DisplacedDijet40_DisplacedTrack_);
    obj.setHT650_DisplacedDijet60_Inclusive(HLT_HT650_DisplacedDijet60_Inclusive_);
    obj.setHT550_DisplacedDijet60_Inclusive(HLT_HT550_DisplacedDijet60_Inclusive_);
    obj.setDiJet110_35_Mjj650_PFMET110(HLT_DiJet110_35_Mjj650_PFMET110_);
    obj.setDiJet110_35_Mjj650_PFMET120(HLT_DiJet110_35_Mjj650_PFMET120_);
    obj.setDiJet110_35_Mjj650_PFMET130(HLT_DiJet110_35_Mjj650_PFMET130_);
    obj.setTripleJet110_35_35_Mjj650_PFMET110(HLT_TripleJet110_35_35_Mjj650_PFMET110_);
    obj.setTripleJet110_35_35_Mjj650_PFMET120(HLT_TripleJet110_35_35_Mjj650_PFMET120_);
    obj.setTripleJet110_35_35_Mjj650_PFMET130(HLT_TripleJet110_35_35_Mjj650_PFMET130_);
    obj.setEle30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned(HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_);
    obj.setEle28_eta2p1_WPTight_Gsf_HT150(HLT_Ele28_eta2p1_WPTight_Gsf_HT150_);
    obj.setEle28_HighEta_SC20_Mass55(HLT_Ele28_HighEta_SC20_Mass55_);
    obj.setDoubleMu20_7_Mass0to30_Photon23(HLT_DoubleMu20_7_Mass0to30_Photon23_);
    obj.setEle15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5(HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_);
    obj.setEle15_IsoVVVL_PFHT450_PFMET50(HLT_Ele15_IsoVVVL_PFHT450_PFMET50_);
    obj.setEle15_IsoVVVL_PFHT450(HLT_Ele15_IsoVVVL_PFHT450_);
    obj.setEle50_IsoVVVL_PFHT450(HLT_Ele50_IsoVVVL_PFHT450_);
    obj.setEle15_IsoVVVL_PFHT600(HLT_Ele15_IsoVVVL_PFHT600_);
    obj.setMu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60(HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_);
    obj.setMu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60(HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_);
    obj.setMu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60(HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_);
    obj.setMu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5(HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_);
    obj.setMu15_IsoVVVL_PFHT450_PFMET50(HLT_Mu15_IsoVVVL_PFHT450_PFMET50_);
    obj.setMu15_IsoVVVL_PFHT450(HLT_Mu15_IsoVVVL_PFHT450_);
    obj.setMu50_IsoVVVL_PFHT450(HLT_Mu50_IsoVVVL_PFHT450_);
    obj.setMu15_IsoVVVL_PFHT600(HLT_Mu15_IsoVVVL_PFHT600_);
    obj.setMu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight(HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_);
    obj.setMu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight(HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_);
    obj.setMu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight(HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_);
    obj.setMu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight(HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_);
    obj.setMu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight(HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_);
    obj.setMu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight(HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_);
    obj.setMu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight(HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_);
    obj.setMu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight(HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_);
    obj.setDimuon10_PsiPrime_Barrel_Seagulls(HLT_Dimuon10_PsiPrime_Barrel_Seagulls_);
    obj.setDimuon20_Jpsi_Barrel_Seagulls(HLT_Dimuon20_Jpsi_Barrel_Seagulls_);
    obj.setDimuon12_Upsilon_y1p4(HLT_Dimuon12_Upsilon_y1p4_);
    obj.setDimuon14_Phi_Barrel_Seagulls(HLT_Dimuon14_Phi_Barrel_Seagulls_);
    obj.setDimuon18_PsiPrime(HLT_Dimuon18_PsiPrime_);
    obj.setDimuon25_Jpsi(HLT_Dimuon25_Jpsi_);
    obj.setDimuon18_PsiPrime_noCorrL1(HLT_Dimuon18_PsiPrime_noCorrL1_);
    obj.setDimuon24_Upsilon_noCorrL1(HLT_Dimuon24_Upsilon_noCorrL1_);
    obj.setDimuon24_Phi_noCorrL1(HLT_Dimuon24_Phi_noCorrL1_);
    obj.setDimuon25_Jpsi_noCorrL1(HLT_Dimuon25_Jpsi_noCorrL1_);
    obj.setDiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8(HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_);
    obj.setDiMu9_Ele9_CaloIdL_TrackIdL_DZ(HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_);
    obj.setDiMu9_Ele9_CaloIdL_TrackIdL(HLT_DiMu9_Ele9_CaloIdL_TrackIdL_);
    obj.setDoubleIsoMu20_eta2p1(HLT_DoubleIsoMu20_eta2p1_);
    obj.setTrkMu12_DoubleTrkMu5NoFiltersNoVtx(HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_);
    obj.setTrkMu16_DoubleTrkMu6NoFiltersNoVtx(HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_);
    obj.setTrkMu17_DoubleTrkMu8NoFiltersNoVtx(HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_);
    obj.setMu8(HLT_Mu8_);
    obj.setMu17(HLT_Mu17_);
    obj.setMu19(HLT_Mu19_);
    obj.setMu17_Photon30_IsoCaloId(HLT_Mu17_Photon30_IsoCaloId_);
    obj.setEle8_CaloIdL_TrackIdL_IsoVL_PFJet30(HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_);
    obj.setEle12_CaloIdL_TrackIdL_IsoVL_PFJet30(HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_);
    obj.setEle15_CaloIdL_TrackIdL_IsoVL_PFJet30(HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_);
    obj.setEle23_CaloIdL_TrackIdL_IsoVL_PFJet30(HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_);
    obj.setEle8_CaloIdM_TrackIdM_PFJet30(HLT_Ele8_CaloIdM_TrackIdM_PFJet30_);
    obj.setEle17_CaloIdM_TrackIdM_PFJet30(HLT_Ele17_CaloIdM_TrackIdM_PFJet30_);
    obj.setEle23_CaloIdM_TrackIdM_PFJet30(HLT_Ele23_CaloIdM_TrackIdM_PFJet30_);
    obj.setEle50_CaloIdVT_GsfTrkIdT_PFJet165(HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_);
    obj.setEle115_CaloIdVT_GsfTrkIdT(HLT_Ele115_CaloIdVT_GsfTrkIdT_);
    obj.setEle135_CaloIdVT_GsfTrkIdT(HLT_Ele135_CaloIdVT_GsfTrkIdT_);
    obj.setEle145_CaloIdVT_GsfTrkIdT(HLT_Ele145_CaloIdVT_GsfTrkIdT_);
    obj.setEle200_CaloIdVT_GsfTrkIdT(HLT_Ele200_CaloIdVT_GsfTrkIdT_);
    obj.setEle250_CaloIdVT_GsfTrkIdT(HLT_Ele250_CaloIdVT_GsfTrkIdT_);
    obj.setEle300_CaloIdVT_GsfTrkIdT(HLT_Ele300_CaloIdVT_GsfTrkIdT_);
    obj.setPFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5(HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_);
    obj.setPFHT330PT30_QuadPFJet_75_60_45_40(HLT_PFHT330PT30_QuadPFJet_75_60_45_40_);
    obj.setPFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94(HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_);
    obj.setPFHT400_SixPFJet32(HLT_PFHT400_SixPFJet32_);
    obj.setPFHT450_SixPFJet36_PFBTagDeepCSV_1p59(HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_);
    obj.setPFHT450_SixPFJet36(HLT_PFHT450_SixPFJet36_);
    obj.setPFHT350(HLT_PFHT350_);
    obj.setPFHT350MinPFJet15(HLT_PFHT350MinPFJet15_);
    obj.setPhoton60_R9Id90_CaloIdL_IsoL(HLT_Photon60_R9Id90_CaloIdL_IsoL_);
    obj.setPhoton60_R9Id90_CaloIdL_IsoL_DisplacedIdL(HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_);
    obj.setPhoton60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15(HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_);
    obj.setECALHT800(HLT_ECALHT800_);
    obj.setDiSC30_18_EIso_AND_HE_Mass70(HLT_DiSC30_18_EIso_AND_HE_Mass70_);
    obj.setPhysics(HLT_Physics_);
    obj.setPhysics_part0(HLT_Physics_part0_);
    obj.setPhysics_part1(HLT_Physics_part1_);
    obj.setPhysics_part2(HLT_Physics_part2_);
    obj.setPhysics_part3(HLT_Physics_part3_);
    obj.setPhysics_part4(HLT_Physics_part4_);
    obj.setPhysics_part5(HLT_Physics_part5_);
    obj.setPhysics_part6(HLT_Physics_part6_);
    obj.setPhysics_part7(HLT_Physics_part7_);
    obj.setRandom(HLT_Random_);
    obj.setZeroBias(HLT_ZeroBias_);
    obj.setZeroBias_Alignment(HLT_ZeroBias_Alignment_);
    obj.setZeroBias_part0(HLT_ZeroBias_part0_);
    obj.setZeroBias_part1(HLT_ZeroBias_part1_);
    obj.setZeroBias_part2(HLT_ZeroBias_part2_);
    obj.setZeroBias_part3(HLT_ZeroBias_part3_);
    obj.setZeroBias_part4(HLT_ZeroBias_part4_);
    obj.setZeroBias_part5(HLT_ZeroBias_part5_);
    obj.setZeroBias_part6(HLT_ZeroBias_part6_);
    obj.setZeroBias_part7(HLT_ZeroBias_part7_);
    obj.setAK4CaloJet30(HLT_AK4CaloJet30_);
    obj.setAK4CaloJet40(HLT_AK4CaloJet40_);
    obj.setAK4CaloJet50(HLT_AK4CaloJet50_);
    obj.setAK4CaloJet80(HLT_AK4CaloJet80_);
    obj.setAK4CaloJet100(HLT_AK4CaloJet100_);
    obj.setAK4CaloJet120(HLT_AK4CaloJet120_);
    obj.setAK4PFJet30(HLT_AK4PFJet30_);
    obj.setAK4PFJet50(HLT_AK4PFJet50_);
    obj.setAK4PFJet80(HLT_AK4PFJet80_);
    obj.setAK4PFJet100(HLT_AK4PFJet100_);
    obj.setAK4PFJet120(HLT_AK4PFJet120_);
    obj.setSinglePhoton10_Eta3p1ForPPRef(HLT_SinglePhoton10_Eta3p1ForPPRef_);
    obj.setSinglePhoton20_Eta3p1ForPPRef(HLT_SinglePhoton20_Eta3p1ForPPRef_);
    obj.setSinglePhoton30_Eta3p1ForPPRef(HLT_SinglePhoton30_Eta3p1ForPPRef_);
    obj.setPhoton20_HoverELoose(HLT_Photon20_HoverELoose_);
    obj.setPhoton30_HoverELoose(HLT_Photon30_HoverELoose_);
    obj.setEcalCalibration(HLT_EcalCalibration_);
    obj.setHcalCalibration(HLT_HcalCalibration_);
    obj.setL1UnpairedBunchBptxMinus(HLT_L1UnpairedBunchBptxMinus_);
    obj.setL1UnpairedBunchBptxPlus(HLT_L1UnpairedBunchBptxPlus_);
    obj.setL1NotBptxOR(HLT_L1NotBptxOR_);
    obj.setL1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142(HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_);
    obj.setCDC_L2cosmic_5_er1p0(HLT_CDC_L2cosmic_5_er1p0_);
    obj.setCDC_L2cosmic_5p5_er1p0(HLT_CDC_L2cosmic_5p5_er1p0_);
    obj.setHcalNZS(HLT_HcalNZS_);
    obj.setHcalPhiSym(HLT_HcalPhiSym_);
    obj.setHcalIsolatedbunch(HLT_HcalIsolatedbunch_);
    obj.setIsoTrackHB(HLT_IsoTrackHB_);
    obj.setIsoTrackHE(HLT_IsoTrackHE_);
    obj.setZeroBias_FirstCollisionAfterAbortGap(HLT_ZeroBias_FirstCollisionAfterAbortGap_);
    obj.setZeroBias_IsolatedBunches(HLT_ZeroBias_IsolatedBunches_);
    obj.setZeroBias_FirstCollisionInTrain(HLT_ZeroBias_FirstCollisionInTrain_);
    obj.setZeroBias_LastCollisionInTrain(HLT_ZeroBias_LastCollisionInTrain_);
    obj.setZeroBias_FirstBXAfterTrain(HLT_ZeroBias_FirstBXAfterTrain_);
    obj.setIsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr(HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_);
    obj.setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_);
    obj.setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_);
    obj.setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_);
    obj.setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_);
    obj.setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_);
    obj.setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_);
    obj.setMediumChargedIsoPFTau50_Trk30_eta2p1_1pr(HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_);
    obj.setMediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr(HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_);
    obj.setMediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1(HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_);
    obj.setMediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1(HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_);
    obj.setMediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1(HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_);
    obj.setEle16_Ele12_Ele8_CaloIdL_TrackIdL(HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_);
    obj.setRsq0p35(HLT_Rsq0p35_);
    obj.setRsq0p40(HLT_Rsq0p40_);
    obj.setRsqMR300_Rsq0p09_MR200(HLT_RsqMR300_Rsq0p09_MR200_);
    obj.setRsqMR320_Rsq0p09_MR200(HLT_RsqMR320_Rsq0p09_MR200_);
    obj.setRsqMR300_Rsq0p09_MR200_4jet(HLT_RsqMR300_Rsq0p09_MR200_4jet_);
    obj.setRsqMR320_Rsq0p09_MR200_4jet(HLT_RsqMR320_Rsq0p09_MR200_4jet_);
    obj.setIsoMu27_MET90(HLT_IsoMu27_MET90_);
    obj.setDoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg(HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_);
    obj.setDoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg(HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_);
    obj.setDoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg(HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_);
    obj.setDoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg(HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_);
    obj.setDoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg(HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_);
    obj.setDoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg(HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_);
    obj.setDoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg(HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_);
    obj.setDoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg(HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_);
    obj.setVBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1(HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_);
    obj.setVBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1(HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_);
    obj.setVBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1(HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_);
    obj.setPhoton50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50(HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_);
    obj.setPhoton75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3(HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_);
    obj.setPhoton75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3(HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_);
    obj.setPFMET100_PFMHT100_IDTight_PFHT60(HLT_PFMET100_PFMHT100_IDTight_PFHT60_);
    obj.setPFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60(HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_);
    obj.setPFMETTypeOne100_PFMHT100_IDTight_PFHT60(HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_);
    obj.setMu18_Mu9_SameSign(HLT_Mu18_Mu9_SameSign_);
    obj.setMu18_Mu9_SameSign_DZ(HLT_Mu18_Mu9_SameSign_DZ_);
    obj.setMu18_Mu9(HLT_Mu18_Mu9_);
    obj.setMu18_Mu9_DZ(HLT_Mu18_Mu9_DZ_);
    obj.setMu20_Mu10_SameSign(HLT_Mu20_Mu10_SameSign_);
    obj.setMu20_Mu10_SameSign_DZ(HLT_Mu20_Mu10_SameSign_DZ_);
    obj.setMu20_Mu10(HLT_Mu20_Mu10_);
    obj.setMu20_Mu10_DZ(HLT_Mu20_Mu10_DZ_);
    obj.setMu23_Mu12_SameSign(HLT_Mu23_Mu12_SameSign_);
    obj.setMu23_Mu12_SameSign_DZ(HLT_Mu23_Mu12_SameSign_DZ_);
    obj.setMu23_Mu12(HLT_Mu23_Mu12_);
    obj.setMu23_Mu12_DZ(HLT_Mu23_Mu12_DZ_);
    obj.setDoubleMu2_Jpsi_DoubleTrk1_Phi1p05(HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_);
    obj.setDoubleMu2_Jpsi_DoubleTkMu0_Phi(HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_);
    obj.setDoubleMu3_DCA_PFMET50_PFMHT60(HLT_DoubleMu3_DCA_PFMET50_PFMHT60_);
    obj.setTripleMu_5_3_3_Mass3p8_DCA(HLT_TripleMu_5_3_3_Mass3p8_DCA_);
    obj.setQuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1(HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_);
    obj.setQuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1(HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_);
    obj.setQuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1(HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_);
    obj.setQuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2(HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_);
    obj.setQuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2(HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_);
    obj.setQuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2(HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_);
    obj.setQuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2(HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_);
    obj.setQuadPFJet98_83_71_15(HLT_QuadPFJet98_83_71_15_);
    obj.setQuadPFJet103_88_75_15(HLT_QuadPFJet103_88_75_15_);
    obj.setQuadPFJet105_88_76_15(HLT_QuadPFJet105_88_76_15_);
    obj.setQuadPFJet111_90_80_15(HLT_QuadPFJet111_90_80_15_);
    obj.setAK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17(HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_);
    obj.setAK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1(HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_);
    obj.setAK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_);
    obj.setAK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_);
    obj.setAK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_);
    obj.setDiphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55(HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_);
    obj.setDiphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto(HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_);
    obj.setMu12_IP6_part0(HLT_Mu12_IP6_part0_);
    obj.setMu12_IP6_part1(HLT_Mu12_IP6_part1_);
    obj.setMu12_IP6_part2(HLT_Mu12_IP6_part2_);
    obj.setMu12_IP6_part3(HLT_Mu12_IP6_part3_);
    obj.setMu12_IP6_part4(HLT_Mu12_IP6_part4_);
    obj.setMu9_IP5_part0(HLT_Mu9_IP5_part0_);
    obj.setMu9_IP5_part1(HLT_Mu9_IP5_part1_);
    obj.setMu9_IP5_part2(HLT_Mu9_IP5_part2_);
    obj.setMu9_IP5_part3(HLT_Mu9_IP5_part3_);
    obj.setMu9_IP5_part4(HLT_Mu9_IP5_part4_);
    obj.setMu7_IP4_part0(HLT_Mu7_IP4_part0_);
    obj.setMu7_IP4_part1(HLT_Mu7_IP4_part1_);
    obj.setMu7_IP4_part2(HLT_Mu7_IP4_part2_);
    obj.setMu7_IP4_part3(HLT_Mu7_IP4_part3_);
    obj.setMu7_IP4_part4(HLT_Mu7_IP4_part4_);
    obj.setMu9_IP4_part0(HLT_Mu9_IP4_part0_);
    obj.setMu9_IP4_part1(HLT_Mu9_IP4_part1_);
    obj.setMu9_IP4_part2(HLT_Mu9_IP4_part2_);
    obj.setMu9_IP4_part3(HLT_Mu9_IP4_part3_);
    obj.setMu9_IP4_part4(HLT_Mu9_IP4_part4_);
    obj.setMu8_IP5_part0(HLT_Mu8_IP5_part0_);
    obj.setMu8_IP5_part1(HLT_Mu8_IP5_part1_);
    obj.setMu8_IP5_part2(HLT_Mu8_IP5_part2_);
    obj.setMu8_IP5_part3(HLT_Mu8_IP5_part3_);
    obj.setMu8_IP5_part4(HLT_Mu8_IP5_part4_);
    obj.setMu8_IP6_part0(HLT_Mu8_IP6_part0_);
    obj.setMu8_IP6_part1(HLT_Mu8_IP6_part1_);
    obj.setMu8_IP6_part2(HLT_Mu8_IP6_part2_);
    obj.setMu8_IP6_part3(HLT_Mu8_IP6_part3_);
    obj.setMu8_IP6_part4(HLT_Mu8_IP6_part4_);
    obj.setMu9_IP6_part0(HLT_Mu9_IP6_part0_);
    obj.setMu9_IP6_part1(HLT_Mu9_IP6_part1_);
    obj.setMu9_IP6_part2(HLT_Mu9_IP6_part2_);
    obj.setMu9_IP6_part3(HLT_Mu9_IP6_part3_);
    obj.setMu9_IP6_part4(HLT_Mu9_IP6_part4_);
    obj.setMu8_IP3_part0(HLT_Mu8_IP3_part0_);
    obj.setMu8_IP3_part1(HLT_Mu8_IP3_part1_);
    obj.setMu8_IP3_part2(HLT_Mu8_IP3_part2_);
    obj.setMu8_IP3_part3(HLT_Mu8_IP3_part3_);
    obj.setMu8_IP3_part4(HLT_Mu8_IP3_part4_);
    obj.setQuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1(HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_);
    obj.setTrkMu6NoFiltersNoVtx(HLT_TrkMu6NoFiltersNoVtx_);
    obj.setTrkMu16NoFiltersNoVtx(HLT_TrkMu16NoFiltersNoVtx_);
    obj.setDoubleTrkMu_16_6_NoFiltersNoVtx(HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_);
    obj.setHLTriggerFinalPath(HLTriggerFinalPath);
    
  
    return obj;
  }
  
  const vector<Svs>& svs(){
    if( SV_.size() > 0) return SV_;
    loadSvs();
  	SV_.reserve(nSV);
    for(size_t idx = 0; idx < nSV; ++idx ){
      Svs obj;
      obj.setdlen(SV_dlen_[idx]);
      obj.setdlenSig(SV_dlenSig_[idx]);
      obj.setpAngle(SV_pAngle_[idx]);
      obj.setchi2(SV_chi2_[idx]);
      obj.setndof(SV_ndof_[idx]);
      obj.setx(SV_x_[idx]);
      obj.sety(SV_y_[idx]);
      obj.setz(SV_z_[idx]);
      obj.setLorentzVector(SV_pt_[idx], SV_eta_[idx], SV_phi_[idx], SV_mass_[idx]);
      SV_.push_back( obj );
    }
    return SV_;
  }
  
  const vector<Electrons>& electrons(){
    if( Electron_.size() > 0) return Electron_;
    loadElectrons();
  	Electron_.reserve(nElectron);
    for(size_t idx = 0; idx < nElectron; ++idx ){
      Electrons obj;
      obj.setdeltaEtaSC(Electron_deltaEtaSC_[idx]);
      obj.setdr03EcalRecHitSumEt(Electron_dr03EcalRecHitSumEt_[idx]);
      obj.setdr03HcalDepth1TowerSumEt(Electron_dr03HcalDepth1TowerSumEt_[idx]);
      obj.setdr03TkSumPt(Electron_dr03TkSumPt_[idx]);
      obj.setdxy(Electron_dxy_[idx]);
      obj.setdxyErr(Electron_dxyErr_[idx]);
      obj.setdz(Electron_dz_[idx]);
      obj.setdzErr(Electron_dzErr_[idx]);
      obj.seteInvMinusPInv(Electron_eInvMinusPInv_[idx]);
      obj.setenergyErr(Electron_energyErr_[idx]);
      obj.sethoe(Electron_hoe_[idx]);
      obj.setip3d(Electron_ip3d_[idx]);
      obj.setminiPFRelIso_all(Electron_miniPFRelIso_all_[idx]);
      obj.setminiPFRelIso_chg(Electron_miniPFRelIso_chg_[idx]);
      obj.setmvaFall17V1Iso(Electron_mvaFall17V1Iso_[idx]);
      obj.setmvaFall17V1noIso(Electron_mvaFall17V1noIso_[idx]);
      obj.setmvaFall17V2Iso(Electron_mvaFall17V2Iso_[idx]);
      obj.setmvaFall17V2noIso(Electron_mvaFall17V2noIso_[idx]);
      obj.setpfRelIso03_all(Electron_pfRelIso03_all_[idx]);
      obj.setpfRelIso03_chg(Electron_pfRelIso03_chg_[idx]);
      obj.setr9(Electron_r9_[idx]);
      obj.setsieie(Electron_sieie_[idx]);
      obj.setsip3d(Electron_sip3d_[idx]);
      obj.setmvaTTH(Electron_mvaTTH_[idx]);
      obj.setcharge(Electron_charge_[idx]);
      obj.setcutBased(Electron_cutBased_[idx]);
      obj.setcutBased_Fall17_V1(Electron_cutBased_Fall17_V1_[idx]);
      obj.setjetIdx(Electron_jetIdx_[idx]);
      obj.setpdgId(Electron_pdgId_[idx]);
      obj.setphotonIdx(Electron_photonIdx_[idx]);
      obj.settightCharge(Electron_tightCharge_[idx]);
      obj.setvidNestedWPBitmap(Electron_vidNestedWPBitmap_[idx]);
      obj.setconvVeto(Electron_convVeto_[idx]);
      obj.setcutBased_HEEP(Electron_cutBased_HEEP_[idx]);
      obj.setisPFcand(Electron_isPFcand_[idx]);
      obj.setlostHits(Electron_lostHits_[idx]);
      obj.setmvaFall17V1Iso_WP80(Electron_mvaFall17V1Iso_WP80_[idx]);
      obj.setmvaFall17V1Iso_WP90(Electron_mvaFall17V1Iso_WP90_[idx]);
      obj.setmvaFall17V1Iso_WPL(Electron_mvaFall17V1Iso_WPL_[idx]);
      obj.setmvaFall17V1noIso_WP80(Electron_mvaFall17V1noIso_WP80_[idx]);
      obj.setmvaFall17V1noIso_WP90(Electron_mvaFall17V1noIso_WP90_[idx]);
      obj.setmvaFall17V1noIso_WPL(Electron_mvaFall17V1noIso_WPL_[idx]);
      obj.setmvaFall17V2Iso_WP80(Electron_mvaFall17V2Iso_WP80_[idx]);
      obj.setmvaFall17V2Iso_WP90(Electron_mvaFall17V2Iso_WP90_[idx]);
      obj.setmvaFall17V2Iso_WPL(Electron_mvaFall17V2Iso_WPL_[idx]);
      obj.setmvaFall17V2noIso_WP80(Electron_mvaFall17V2noIso_WP80_[idx]);
      obj.setmvaFall17V2noIso_WP90(Electron_mvaFall17V2noIso_WP90_[idx]);
      obj.setmvaFall17V2noIso_WPL(Electron_mvaFall17V2noIso_WPL_[idx]);
      obj.setgenPartIdx(Electron_genPartIdx_[idx]);
      obj.setgenPartFlav(Electron_genPartFlav_[idx]);
      obj.setcleanmask(Electron_cleanmask_[idx]);
      obj.setLorentzVector(Electron_pt_[idx], Electron_eta_[idx], Electron_phi_[idx], Electron_mass_[idx]);
      Electron_.push_back( obj );
    }
    return Electron_;
  }
  
  const Met met(){
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
  
  const Genmet genmet(){
    //non-vectorial objects are recomputed every
    //time for simplicity 
    loadGenmet();
  
    Genmet obj;
    
    
  
    return obj;
  }
  
  const vector<Fatjets>& fatjets(){
    if( FatJet_.size() > 0) return FatJet_;
    loadFatjets();
  	FatJet_.reserve(nFatJet);
    for(size_t idx = 0; idx < nFatJet; ++idx ){
      Fatjets obj;
      obj.setarea(FatJet_area_[idx]);
      obj.setbtagCMVA(FatJet_btagCMVA_[idx]);
      obj.setbtagCSVV2(FatJet_btagCSVV2_[idx]);
      obj.setbtagDeepB(FatJet_btagDeepB_[idx]);
      obj.setbtagHbb(FatJet_btagHbb_[idx]);
      obj.setmsoftdrop(FatJet_msoftdrop_[idx]);
      obj.setn2b1(FatJet_n2b1_[idx]);
      obj.setn3b1(FatJet_n3b1_[idx]);
      obj.setrawFactor(FatJet_rawFactor_[idx]);
      obj.settau1(FatJet_tau1_[idx]);
      obj.settau2(FatJet_tau2_[idx]);
      obj.settau3(FatJet_tau3_[idx]);
      obj.settau4(FatJet_tau4_[idx]);
      obj.setjetId(FatJet_jetId_[idx]);
      obj.setsubJetIdx1(FatJet_subJetIdx1_[idx]);
      obj.setsubJetIdx2(FatJet_subJetIdx2_[idx]);
      obj.setLorentzVector(FatJet_pt_[idx], FatJet_eta_[idx], FatJet_phi_[idx], FatJet_mass_[idx]);
      FatJet_.push_back( obj );
    }
    return FatJet_;
  }
  
  const vector<Subjets>& subjets(){
    if( SubJet_.size() > 0) return SubJet_;
    loadSubjets();
  	SubJet_.reserve(nSubJet);
    for(size_t idx = 0; idx < nSubJet; ++idx ){
      Subjets obj;
      obj.setbtagCMVA(SubJet_btagCMVA_[idx]);
      obj.setbtagCSVV2(SubJet_btagCSVV2_[idx]);
      obj.setbtagDeepB(SubJet_btagDeepB_[idx]);
      obj.setn2b1(SubJet_n2b1_[idx]);
      obj.setn3b1(SubJet_n3b1_[idx]);
      obj.settau1(SubJet_tau1_[idx]);
      obj.settau2(SubJet_tau2_[idx]);
      obj.settau3(SubJet_tau3_[idx]);
      obj.settau4(SubJet_tau4_[idx]);
      obj.setLorentzVector(SubJet_pt_[idx], SubJet_eta_[idx], SubJet_phi_[idx], SubJet_mass_[idx]);
      SubJet_.push_back( obj );
    }
    return SubJet_;
  }
  
  const vector<Subgenjetak8s>& subgenjetak8s(){
    if( SubGenJetAK8_.size() > 0) return SubGenJetAK8_;
    loadSubgenjetak8s();
  	SubGenJetAK8_.reserve(nSubGenJetAK8);
    for(size_t idx = 0; idx < nSubGenJetAK8; ++idx ){
      Subgenjetak8s obj;
      
      obj.setLorentzVector(SubGenJetAK8_pt_[idx], SubGenJetAK8_eta_[idx], SubGenJetAK8_phi_[idx], SubGenJetAK8_mass_[idx]);
      SubGenJetAK8_.push_back( obj );
    }
    return SubGenJetAK8_;
  }
  
  const Pileup pileup(){
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
  
  const Lheweight lheweight(){
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
  Float_t btagWeight_CSVV2_;
  Float_t btagWeight_DeepCSVB_;
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
  Float_t *Electron_eInvMinusPInv_;
  Float_t *Electron_energyErr_;
  Float_t *Electron_eta_;
  Float_t *Electron_hoe_;
  Float_t *Electron_ip3d_;
  Float_t *Electron_mass_;
  Float_t *Electron_miniPFRelIso_all_;
  Float_t *Electron_miniPFRelIso_chg_;
  Float_t *Electron_mvaFall17V1Iso_;
  Float_t *Electron_mvaFall17V1noIso_;
  Float_t *Electron_mvaFall17V2Iso_;
  Float_t *Electron_mvaFall17V2noIso_;
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
  Int_t *Electron_cutBased_Fall17_V1_;
  Int_t *Electron_jetIdx_;
  Int_t *Electron_pdgId_;
  Int_t *Electron_photonIdx_;
  Int_t *Electron_tightCharge_;
  Int_t *Electron_vidNestedWPBitmap_;
  Bool_t *Electron_convVeto_;
  Bool_t *Electron_cutBased_HEEP_;
  Bool_t *Electron_isPFcand_;
  UChar_t *Electron_lostHits_;
  Bool_t *Electron_mvaFall17V1Iso_WP80_;
  Bool_t *Electron_mvaFall17V1Iso_WP90_;
  Bool_t *Electron_mvaFall17V1Iso_WPL_;
  Bool_t *Electron_mvaFall17V1noIso_WP80_;
  Bool_t *Electron_mvaFall17V1noIso_WP90_;
  Bool_t *Electron_mvaFall17V1noIso_WPL_;
  Bool_t *Electron_mvaFall17V2Iso_WP80_;
  Bool_t *Electron_mvaFall17V2Iso_WP90_;
  Bool_t *Electron_mvaFall17V2Iso_WPL_;
  Bool_t *Electron_mvaFall17V2noIso_WP80_;
  Bool_t *Electron_mvaFall17V2noIso_WP90_;
  Bool_t *Electron_mvaFall17V2noIso_WPL_;
  Float_t *FatJet_area_;
  Float_t *FatJet_btagCMVA_;
  Float_t *FatJet_btagCSVV2_;
  Float_t *FatJet_btagDeepB_;
  Float_t *FatJet_btagHbb_;
  Float_t *FatJet_eta_;
  Float_t *FatJet_mass_;
  Float_t *FatJet_msoftdrop_;
  Float_t *FatJet_n2b1_;
  Float_t *FatJet_n3b1_;
  Float_t *FatJet_phi_;
  Float_t *FatJet_pt_;
  Float_t *FatJet_rawFactor_;
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
  Float_t *SubGenJetAK8_eta_;
  Float_t *SubGenJetAK8_mass_;
  Float_t *SubGenJetAK8_phi_;
  Float_t *SubGenJetAK8_pt_;
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
  Float_t *IsoTrack_dxy_;
  Float_t *IsoTrack_dz_;
  Float_t *IsoTrack_eta_;
  Float_t *IsoTrack_pfRelIso03_all_;
  Float_t *IsoTrack_pfRelIso03_chg_;
  Float_t *IsoTrack_phi_;
  Float_t *IsoTrack_pt_;
  Float_t *IsoTrack_miniPFRelIso_all_;
  Float_t *IsoTrack_miniPFRelIso_chg_;
  Int_t *IsoTrack_pdgId_;
  Bool_t *IsoTrack_isHighPurityTrack_;
  Bool_t *IsoTrack_isPFcand_;
  Float_t *Jet_area_;
  Float_t *Jet_btagCMVA_;
  Float_t *Jet_btagCSVV2_;
  Float_t *Jet_btagDeepB_;
  Float_t *Jet_btagDeepC_;
  Float_t *Jet_btagDeepFlavB_;
  Float_t *Jet_chEmEF_;
  Float_t *Jet_chHEF_;
  Float_t *Jet_eta_;
  Float_t *Jet_mass_;
  Float_t *Jet_muEF_;
  Float_t *Jet_neEmEF_;
  Float_t *Jet_neHEF_;
  Float_t *Jet_phi_;
  Float_t *Jet_pt_;
  Float_t *Jet_qgl_;
  Float_t *Jet_rawFactor_;
  Float_t *Jet_bRegCorr_;
  Float_t *Jet_bRegRes_;
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
  Float_t *LHEPart_pt_;
  Float_t *LHEPart_eta_;
  Float_t *LHEPart_phi_;
  Float_t *LHEPart_mass_;
  Int_t *LHEPart_pdgId_;
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
  Int_t *Photon_cutBasedBitmap_;
  Int_t *Photon_electronIdx_;
  Int_t *Photon_jetIdx_;
  Int_t *Photon_pdgId_;
  Int_t *Photon_vidNestedWPBitmap_;
  Bool_t *Photon_electronVeto_;
  Bool_t *Photon_isScEtaEB_;
  Bool_t *Photon_isScEtaEE_;
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
  Float_t *Tau_rawIsodR03_;
  Int_t *Tau_charge_;
  Int_t *Tau_decayMode_;
  Int_t *Tau_jetIdx_;
  Int_t *Tau_rawAntiEleCat_;
  UChar_t *Tau_idAntiEle_;
  UChar_t *Tau_idAntiMu_;
  Bool_t *Tau_idDecayMode_;
  Bool_t *Tau_idDecayModeNewDMs_;
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
  Bool_t HLT_AK8PFJet380_TrimMass30_;
  Bool_t HLT_AK8PFJet400_TrimMass30_;
  Bool_t HLT_AK8PFJet420_TrimMass30_;
  Bool_t HLT_AK8PFHT750_TrimMass50_;
  Bool_t HLT_AK8PFHT800_TrimMass50_;
  Bool_t HLT_AK8PFHT850_TrimMass50_;
  Bool_t HLT_AK8PFHT900_TrimMass50_;
  Bool_t HLT_CaloJet500_NoJetID_;
  Bool_t HLT_CaloJet550_NoJetID_;
  Bool_t HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_;
  Bool_t HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_;
  Bool_t HLT_Trimuon5_3p5_2_Upsilon_Muon_;
  Bool_t HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_;
  Bool_t HLT_DoubleEle25_CaloIdL_MW_;
  Bool_t HLT_DoubleEle27_CaloIdL_MW_;
  Bool_t HLT_DoubleEle33_CaloIdL_MW_;
  Bool_t HLT_DoubleEle24_eta2p1_WPTight_Gsf_;
  Bool_t HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_;
  Bool_t HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_;
  Bool_t HLT_Ele27_Ele37_CaloIdL_MW_;
  Bool_t HLT_Mu27_Ele37_CaloIdL_MW_;
  Bool_t HLT_Mu37_Ele27_CaloIdL_MW_;
  Bool_t HLT_Mu37_TkMu27_;
  Bool_t HLT_DoubleMu4_3_Bs_;
  Bool_t HLT_DoubleMu4_3_Jpsi_;
  Bool_t HLT_DoubleMu4_JpsiTrk_Displaced_;
  Bool_t HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_;
  Bool_t HLT_DoubleMu3_Trk_Tau3mu_;
  Bool_t HLT_DoubleMu3_TkMu_DsTau3Mu_;
  Bool_t HLT_DoubleMu4_PsiPrimeTrk_Displaced_;
  Bool_t HLT_DoubleMu4_Mass3p8_DZ_PFHT350_;
  Bool_t HLT_Mu3_PFJet40_;
  Bool_t HLT_Mu7p5_L2Mu2_Jpsi_;
  Bool_t HLT_Mu7p5_L2Mu2_Upsilon_;
  Bool_t HLT_Mu7p5_Track2_Jpsi_;
  Bool_t HLT_Mu7p5_Track3p5_Jpsi_;
  Bool_t HLT_Mu7p5_Track7_Jpsi_;
  Bool_t HLT_Mu7p5_Track2_Upsilon_;
  Bool_t HLT_Mu7p5_Track3p5_Upsilon_;
  Bool_t HLT_Mu7p5_Track7_Upsilon_;
  Bool_t HLT_Mu3_L1SingleMu5orSingleMu7_;
  Bool_t HLT_DoublePhoton33_CaloIdL_;
  Bool_t HLT_DoublePhoton70_;
  Bool_t HLT_DoublePhoton85_;
  Bool_t HLT_Ele20_WPTight_Gsf_;
  Bool_t HLT_Ele15_WPLoose_Gsf_;
  Bool_t HLT_Ele17_WPLoose_Gsf_;
  Bool_t HLT_Ele20_WPLoose_Gsf_;
  Bool_t HLT_Ele20_eta2p1_WPLoose_Gsf_;
  Bool_t HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_;
  Bool_t HLT_Ele27_WPTight_Gsf_;
  Bool_t HLT_Ele28_WPTight_Gsf_;
  Bool_t HLT_Ele30_WPTight_Gsf_;
  Bool_t HLT_Ele32_WPTight_Gsf_;
  Bool_t HLT_Ele35_WPTight_Gsf_;
  Bool_t HLT_Ele35_WPTight_Gsf_L1EGMT_;
  Bool_t HLT_Ele38_WPTight_Gsf_;
  Bool_t HLT_Ele40_WPTight_Gsf_;
  Bool_t HLT_Ele32_WPTight_Gsf_L1DoubleEG_;
  Bool_t HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_;
  Bool_t HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_;
  Bool_t HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_;
  Bool_t HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;
  Bool_t HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;
  Bool_t HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_;
  Bool_t HLT_HT450_Beamspot_;
  Bool_t HLT_HT300_Beamspot_;
  Bool_t HLT_ZeroBias_Beamspot_;
  Bool_t HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_;
  Bool_t HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_;
  Bool_t HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_;
  Bool_t HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;
  Bool_t HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;
  Bool_t HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_;
  Bool_t HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_;
  Bool_t HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_;
  Bool_t HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_;
  Bool_t HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_;
  Bool_t HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;
  Bool_t HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;
  Bool_t HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_;
  Bool_t HLT_IsoMu20_;
  Bool_t HLT_IsoMu24_;
  Bool_t HLT_IsoMu24_eta2p1_;
  Bool_t HLT_IsoMu27_;
  Bool_t HLT_IsoMu30_;
  Bool_t HLT_UncorrectedJetE30_NoBPTX_;
  Bool_t HLT_UncorrectedJetE30_NoBPTX3BX_;
  Bool_t HLT_UncorrectedJetE60_NoBPTX3BX_;
  Bool_t HLT_UncorrectedJetE70_NoBPTX3BX_;
  Bool_t HLT_L1SingleMu18_;
  Bool_t HLT_L1SingleMu25_;
  Bool_t HLT_L2Mu10_;
  Bool_t HLT_L2Mu10_NoVertex_NoBPTX3BX_;
  Bool_t HLT_L2Mu10_NoVertex_NoBPTX_;
  Bool_t HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_;
  Bool_t HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_;
  Bool_t HLT_L2Mu50_;
  Bool_t HLT_L2Mu23NoVtx_2Cha_;
  Bool_t HLT_L2Mu23NoVtx_2Cha_CosmicSeed_;
  Bool_t HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_;
  Bool_t HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_;
  Bool_t HLT_DoubleL2Mu50_;
  Bool_t HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_;
  Bool_t HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_;
  Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_;
  Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_;
  Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_;
  Bool_t HLT_DoubleL2Mu23NoVtx_2Cha_;
  Bool_t HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_;
  Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_;
  Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_;
  Bool_t HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_;
  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_;
  Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_;
  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_;
  Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_;
  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_;
  Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_;
  Bool_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_;
  Bool_t HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_;
  Bool_t HLT_Mu25_TkMu0_Onia_;
  Bool_t HLT_Mu30_TkMu0_Psi_;
  Bool_t HLT_Mu30_TkMu0_Upsilon_;
  Bool_t HLT_Mu20_TkMu0_Phi_;
  Bool_t HLT_Mu25_TkMu0_Phi_;
  Bool_t HLT_Mu12_;
  Bool_t HLT_Mu15_;
  Bool_t HLT_Mu20_;
  Bool_t HLT_Mu27_;
  Bool_t HLT_Mu50_;
  Bool_t HLT_Mu55_;
  Bool_t HLT_OldMu100_;
  Bool_t HLT_TkMu100_;
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
  Bool_t HLT_AK8PFJet15_;
  Bool_t HLT_AK8PFJet25_;
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
  Bool_t HLT_AK8PFJet550_;
  Bool_t HLT_PFJet15_;
  Bool_t HLT_PFJet25_;
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
  Bool_t HLT_PFJet550_;
  Bool_t HLT_PFJetFwd15_;
  Bool_t HLT_PFJetFwd25_;
  Bool_t HLT_PFJetFwd40_;
  Bool_t HLT_PFJetFwd60_;
  Bool_t HLT_PFJetFwd80_;
  Bool_t HLT_PFJetFwd140_;
  Bool_t HLT_PFJetFwd200_;
  Bool_t HLT_PFJetFwd260_;
  Bool_t HLT_PFJetFwd320_;
  Bool_t HLT_PFJetFwd400_;
  Bool_t HLT_PFJetFwd450_;
  Bool_t HLT_PFJetFwd500_;
  Bool_t HLT_AK8PFJetFwd15_;
  Bool_t HLT_AK8PFJetFwd25_;
  Bool_t HLT_AK8PFJetFwd40_;
  Bool_t HLT_AK8PFJetFwd60_;
  Bool_t HLT_AK8PFJetFwd80_;
  Bool_t HLT_AK8PFJetFwd140_;
  Bool_t HLT_AK8PFJetFwd200_;
  Bool_t HLT_AK8PFJetFwd260_;
  Bool_t HLT_AK8PFJetFwd320_;
  Bool_t HLT_AK8PFJetFwd400_;
  Bool_t HLT_AK8PFJetFwd450_;
  Bool_t HLT_AK8PFJetFwd500_;
  Bool_t HLT_PFHT180_;
  Bool_t HLT_PFHT250_;
  Bool_t HLT_PFHT370_;
  Bool_t HLT_PFHT430_;
  Bool_t HLT_PFHT510_;
  Bool_t HLT_PFHT590_;
  Bool_t HLT_PFHT680_;
  Bool_t HLT_PFHT780_;
  Bool_t HLT_PFHT890_;
  Bool_t HLT_PFHT1050_;
  Bool_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_;
  Bool_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_;
  Bool_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_;
  Bool_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_;
  Bool_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_;
  Bool_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_;
  Bool_t HLT_PFMET110_PFMHT110_IDTight_;
  Bool_t HLT_PFMET120_PFMHT120_IDTight_;
  Bool_t HLT_PFMET130_PFMHT130_IDTight_;
  Bool_t HLT_PFMET140_PFMHT140_IDTight_;
  Bool_t HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_;
  Bool_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_;
  Bool_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_;
  Bool_t HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_;
  Bool_t HLT_PFMETTypeOne110_PFMHT110_IDTight_;
  Bool_t HLT_PFMETTypeOne120_PFMHT120_IDTight_;
  Bool_t HLT_PFMETTypeOne130_PFMHT130_IDTight_;
  Bool_t HLT_PFMETTypeOne140_PFMHT140_IDTight_;
  Bool_t HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_;
  Bool_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_;
  Bool_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_;
  Bool_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_;
  Bool_t HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_;
  Bool_t HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_;
  Bool_t HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_;
  Bool_t HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_;
  Bool_t HLT_L1ETMHadSeeds_;
  Bool_t HLT_CaloMHT90_;
  Bool_t HLT_CaloMET80_NotCleaned_;
  Bool_t HLT_CaloMET90_NotCleaned_;
  Bool_t HLT_CaloMET100_NotCleaned_;
  Bool_t HLT_CaloMET110_NotCleaned_;
  Bool_t HLT_CaloMET250_NotCleaned_;
  Bool_t HLT_CaloMET70_HBHECleaned_;
  Bool_t HLT_CaloMET80_HBHECleaned_;
  Bool_t HLT_CaloMET90_HBHECleaned_;
  Bool_t HLT_CaloMET100_HBHECleaned_;
  Bool_t HLT_CaloMET250_HBHECleaned_;
  Bool_t HLT_CaloMET300_HBHECleaned_;
  Bool_t HLT_CaloMET350_HBHECleaned_;
  Bool_t HLT_PFMET200_NotCleaned_;
  Bool_t HLT_PFMET200_HBHECleaned_;
  Bool_t HLT_PFMET250_HBHECleaned_;
  Bool_t HLT_PFMET300_HBHECleaned_;
  Bool_t HLT_PFMET200_HBHE_BeamHaloCleaned_;
  Bool_t HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_;
  Bool_t HLT_MET105_IsoTrk50_;
  Bool_t HLT_MET120_IsoTrk50_;
  Bool_t HLT_SingleJet30_Mu12_SinglePFJet40_;
  Bool_t HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_;
  Bool_t HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_;
  Bool_t HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_;
  Bool_t HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_;
  Bool_t HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t HLT_DoublePFJets40_CaloBTagDeepCSV_p71_;
  Bool_t HLT_DoublePFJets100_CaloBTagDeepCSV_p71_;
  Bool_t HLT_DoublePFJets200_CaloBTagDeepCSV_p71_;
  Bool_t HLT_DoublePFJets350_CaloBTagDeepCSV_p71_;
  Bool_t HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_;
  Bool_t HLT_Photon300_NoHE_;
  Bool_t HLT_Mu8_TrkIsoVVL_;
  Bool_t HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_;
  Bool_t HLT_Mu8_DiEle12_CaloIdL_TrackIdL_;
  Bool_t HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_;
  Bool_t HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_;
  Bool_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu17_TrkIsoVVL_;
  Bool_t HLT_Mu19_TrkIsoVVL_;
  Bool_t HLT_BTagMu_AK4DiJet20_Mu5_;
  Bool_t HLT_BTagMu_AK4DiJet40_Mu5_;
  Bool_t HLT_BTagMu_AK4DiJet70_Mu5_;
  Bool_t HLT_BTagMu_AK4DiJet110_Mu5_;
  Bool_t HLT_BTagMu_AK4DiJet170_Mu5_;
  Bool_t HLT_BTagMu_AK4Jet300_Mu5_;
  Bool_t HLT_BTagMu_AK8DiJet170_Mu5_;
  Bool_t HLT_BTagMu_AK8Jet170_DoubleMu5_;
  Bool_t HLT_BTagMu_AK8Jet300_Mu5_;
  Bool_t HLT_BTagMu_AK4DiJet20_Mu5_noalgo_;
  Bool_t HLT_BTagMu_AK4DiJet40_Mu5_noalgo_;
  Bool_t HLT_BTagMu_AK4DiJet70_Mu5_noalgo_;
  Bool_t HLT_BTagMu_AK4DiJet110_Mu5_noalgo_;
  Bool_t HLT_BTagMu_AK4DiJet170_Mu5_noalgo_;
  Bool_t HLT_BTagMu_AK4Jet300_Mu5_noalgo_;
  Bool_t HLT_BTagMu_AK8DiJet170_Mu5_noalgo_;
  Bool_t HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_;
  Bool_t HLT_BTagMu_AK8Jet300_Mu5_noalgo_;
  Bool_t HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_;
  Bool_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_;
  Bool_t HLT_Mu12_DoublePhoton20_;
  Bool_t HLT_TriplePhoton_20_20_20_CaloIdLV2_;
  Bool_t HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_;
  Bool_t HLT_TriplePhoton_30_30_10_CaloIdLV2_;
  Bool_t HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_;
  Bool_t HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_;
  Bool_t HLT_Photon20_;
  Bool_t HLT_Photon33_;
  Bool_t HLT_Photon50_;
  Bool_t HLT_Photon75_;
  Bool_t HLT_Photon90_;
  Bool_t HLT_Photon120_;
  Bool_t HLT_Photon150_;
  Bool_t HLT_Photon175_;
  Bool_t HLT_Photon200_;
  Bool_t HLT_Photon100EB_TightID_TightIso_;
  Bool_t HLT_Photon110EB_TightID_TightIso_;
  Bool_t HLT_Photon120EB_TightID_TightIso_;
  Bool_t HLT_Photon100EBHE10_;
  Bool_t HLT_Photon100EEHE10_;
  Bool_t HLT_Photon100EE_TightID_TightIso_;
  Bool_t HLT_Photon50_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon75_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_;
  Bool_t HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_;
  Bool_t HLT_Photon90_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon120_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon165_R9Id90_HE10_IsoM_;
  Bool_t HLT_Photon90_CaloIdL_PFHT700_;
  Bool_t HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_;
  Bool_t HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_;
  Bool_t HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_;
  Bool_t HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_;
  Bool_t HLT_Photon35_TwoProngs35_;
  Bool_t HLT_IsoMu24_TwoProngs35_;
  Bool_t HLT_Dimuon0_Jpsi_L1_NoOS_;
  Bool_t HLT_Dimuon0_Jpsi_NoVertexing_NoOS_;
  Bool_t HLT_Dimuon0_Jpsi_;
  Bool_t HLT_Dimuon0_Jpsi_NoVertexing_;
  Bool_t HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_;
  Bool_t HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_;
  Bool_t HLT_Dimuon0_Jpsi3p5_Muon2_;
  Bool_t HLT_Dimuon0_Upsilon_L1_4p5_;
  Bool_t HLT_Dimuon0_Upsilon_L1_5_;
  Bool_t HLT_Dimuon0_Upsilon_L1_4p5NoOS_;
  Bool_t HLT_Dimuon0_Upsilon_L1_4p5er2p0_;
  Bool_t HLT_Dimuon0_Upsilon_L1_4p5er2p0M_;
  Bool_t HLT_Dimuon0_Upsilon_NoVertexing_;
  Bool_t HLT_Dimuon0_Upsilon_L1_5M_;
  Bool_t HLT_Dimuon0_LowMass_L1_0er1p5R_;
  Bool_t HLT_Dimuon0_LowMass_L1_0er1p5_;
  Bool_t HLT_Dimuon0_LowMass_;
  Bool_t HLT_Dimuon0_LowMass_L1_4_;
  Bool_t HLT_Dimuon0_LowMass_L1_4R_;
  Bool_t HLT_Dimuon0_LowMass_L1_TM530_;
  Bool_t HLT_Dimuon0_Upsilon_Muon_L1_TM0_;
  Bool_t HLT_Dimuon0_Upsilon_Muon_NoL1Mass_;
  Bool_t HLT_TripleMu_5_3_3_Mass3p8_DZ_;
  Bool_t HLT_TripleMu_10_5_5_DZ_;
  Bool_t HLT_TripleMu_12_10_5_;
  Bool_t HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_;
  Bool_t HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_;
  Bool_t HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_;
  Bool_t HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_;
  Bool_t HLT_DoubleMu3_DZ_PFMET50_PFMHT60_;
  Bool_t HLT_DoubleMu3_DZ_PFMET70_PFMHT70_;
  Bool_t HLT_DoubleMu3_DZ_PFMET90_PFMHT90_;
  Bool_t HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_;
  Bool_t HLT_DoubleMu4_Jpsi_Displaced_;
  Bool_t HLT_DoubleMu4_Jpsi_NoVertexing_;
  Bool_t HLT_DoubleMu4_JpsiTrkTrk_Displaced_;
  Bool_t HLT_DoubleMu43NoFiltersNoVtx_;
  Bool_t HLT_DoubleMu48NoFiltersNoVtx_;
  Bool_t HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_;
  Bool_t HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_;
  Bool_t HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_;
  Bool_t HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_;
  Bool_t HLT_DoubleMu33NoFiltersNoVtxDisplaced_;
  Bool_t HLT_DoubleMu40NoFiltersNoVtxDisplaced_;
  Bool_t HLT_DoubleMu20_7_Mass0to30_L1_DM4_;
  Bool_t HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_;
  Bool_t HLT_HT425_;
  Bool_t HLT_HT430_DisplacedDijet40_DisplacedTrack_;
  Bool_t HLT_HT500_DisplacedDijet40_DisplacedTrack_;
  Bool_t HLT_HT430_DisplacedDijet60_DisplacedTrack_;
  Bool_t HLT_HT400_DisplacedDijet40_DisplacedTrack_;
  Bool_t HLT_HT650_DisplacedDijet60_Inclusive_;
  Bool_t HLT_HT550_DisplacedDijet60_Inclusive_;
  Bool_t HLT_DiJet110_35_Mjj650_PFMET110_;
  Bool_t HLT_DiJet110_35_Mjj650_PFMET120_;
  Bool_t HLT_DiJet110_35_Mjj650_PFMET130_;
  Bool_t HLT_TripleJet110_35_35_Mjj650_PFMET110_;
  Bool_t HLT_TripleJet110_35_35_Mjj650_PFMET120_;
  Bool_t HLT_TripleJet110_35_35_Mjj650_PFMET130_;
  Bool_t HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_;
  Bool_t HLT_Ele28_eta2p1_WPTight_Gsf_HT150_;
  Bool_t HLT_Ele28_HighEta_SC20_Mass55_;
  Bool_t HLT_DoubleMu20_7_Mass0to30_Photon23_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT450_PFMET50_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT450_;
  Bool_t HLT_Ele50_IsoVVVL_PFHT450_;
  Bool_t HLT_Ele15_IsoVVVL_PFHT600_;
  Bool_t HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_;
  Bool_t HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_;
  Bool_t HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT450_PFMET50_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT450_;
  Bool_t HLT_Mu50_IsoVVVL_PFHT450_;
  Bool_t HLT_Mu15_IsoVVVL_PFHT600_;
  Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_;
  Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_;
  Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_;
  Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_;
  Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_;
  Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_;
  Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_;
  Bool_t HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_;
  Bool_t HLT_Dimuon10_PsiPrime_Barrel_Seagulls_;
  Bool_t HLT_Dimuon20_Jpsi_Barrel_Seagulls_;
  Bool_t HLT_Dimuon12_Upsilon_y1p4_;
  Bool_t HLT_Dimuon14_Phi_Barrel_Seagulls_;
  Bool_t HLT_Dimuon18_PsiPrime_;
  Bool_t HLT_Dimuon25_Jpsi_;
  Bool_t HLT_Dimuon18_PsiPrime_noCorrL1_;
  Bool_t HLT_Dimuon24_Upsilon_noCorrL1_;
  Bool_t HLT_Dimuon24_Phi_noCorrL1_;
  Bool_t HLT_Dimuon25_Jpsi_noCorrL1_;
  Bool_t HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8_;
  Bool_t HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_;
  Bool_t HLT_DiMu9_Ele9_CaloIdL_TrackIdL_;
  Bool_t HLT_DoubleIsoMu20_eta2p1_;
  Bool_t HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_;
  Bool_t HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_;
  Bool_t HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_;
  Bool_t HLT_Mu8_;
  Bool_t HLT_Mu17_;
  Bool_t HLT_Mu19_;
  Bool_t HLT_Mu17_Photon30_IsoCaloId_;
  Bool_t HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_;
  Bool_t HLT_Ele8_CaloIdM_TrackIdM_PFJet30_;
  Bool_t HLT_Ele17_CaloIdM_TrackIdM_PFJet30_;
  Bool_t HLT_Ele23_CaloIdM_TrackIdM_PFJet30_;
  Bool_t HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_;
  Bool_t HLT_Ele115_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_Ele135_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_Ele145_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_Ele200_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_Ele250_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_Ele300_CaloIdVT_GsfTrkIdT_;
  Bool_t HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_;
  Bool_t HLT_PFHT330PT30_QuadPFJet_75_60_45_40_;
  Bool_t HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_;
  Bool_t HLT_PFHT400_SixPFJet32_;
  Bool_t HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_;
  Bool_t HLT_PFHT450_SixPFJet36_;
  Bool_t HLT_PFHT350_;
  Bool_t HLT_PFHT350MinPFJet15_;
  Bool_t HLT_Photon60_R9Id90_CaloIdL_IsoL_;
  Bool_t HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_;
  Bool_t HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_;
  Bool_t HLT_ECALHT800_;
  Bool_t HLT_DiSC30_18_EIso_AND_HE_Mass70_;
  Bool_t HLT_Physics_;
  Bool_t HLT_Physics_part0_;
  Bool_t HLT_Physics_part1_;
  Bool_t HLT_Physics_part2_;
  Bool_t HLT_Physics_part3_;
  Bool_t HLT_Physics_part4_;
  Bool_t HLT_Physics_part5_;
  Bool_t HLT_Physics_part6_;
  Bool_t HLT_Physics_part7_;
  Bool_t HLT_Random_;
  Bool_t HLT_ZeroBias_;
  Bool_t HLT_ZeroBias_Alignment_;
  Bool_t HLT_ZeroBias_part0_;
  Bool_t HLT_ZeroBias_part1_;
  Bool_t HLT_ZeroBias_part2_;
  Bool_t HLT_ZeroBias_part3_;
  Bool_t HLT_ZeroBias_part4_;
  Bool_t HLT_ZeroBias_part5_;
  Bool_t HLT_ZeroBias_part6_;
  Bool_t HLT_ZeroBias_part7_;
  Bool_t HLT_AK4CaloJet30_;
  Bool_t HLT_AK4CaloJet40_;
  Bool_t HLT_AK4CaloJet50_;
  Bool_t HLT_AK4CaloJet80_;
  Bool_t HLT_AK4CaloJet100_;
  Bool_t HLT_AK4CaloJet120_;
  Bool_t HLT_AK4PFJet30_;
  Bool_t HLT_AK4PFJet50_;
  Bool_t HLT_AK4PFJet80_;
  Bool_t HLT_AK4PFJet100_;
  Bool_t HLT_AK4PFJet120_;
  Bool_t HLT_SinglePhoton10_Eta3p1ForPPRef_;
  Bool_t HLT_SinglePhoton20_Eta3p1ForPPRef_;
  Bool_t HLT_SinglePhoton30_Eta3p1ForPPRef_;
  Bool_t HLT_Photon20_HoverELoose_;
  Bool_t HLT_Photon30_HoverELoose_;
  Bool_t HLT_EcalCalibration_;
  Bool_t HLT_HcalCalibration_;
  Bool_t HLT_L1UnpairedBunchBptxMinus_;
  Bool_t HLT_L1UnpairedBunchBptxPlus_;
  Bool_t HLT_L1NotBptxOR_;
  Bool_t HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_;
  Bool_t HLT_CDC_L2cosmic_5_er1p0_;
  Bool_t HLT_CDC_L2cosmic_5p5_er1p0_;
  Bool_t HLT_HcalNZS_;
  Bool_t HLT_HcalPhiSym_;
  Bool_t HLT_HcalIsolatedbunch_;
  Bool_t HLT_IsoTrackHB_;
  Bool_t HLT_IsoTrackHE_;
  Bool_t HLT_ZeroBias_FirstCollisionAfterAbortGap_;
  Bool_t HLT_ZeroBias_IsolatedBunches_;
  Bool_t HLT_ZeroBias_FirstCollisionInTrain_;
  Bool_t HLT_ZeroBias_LastCollisionInTrain_;
  Bool_t HLT_ZeroBias_FirstBXAfterTrain_;
  Bool_t HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_;
  Bool_t HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_;
  Bool_t HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_;
  Bool_t HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_;
  Bool_t HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_;
  Bool_t HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_;
  Bool_t HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_;
  Bool_t HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_;
  Bool_t HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_;
  Bool_t HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_;
  Bool_t HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_;
  Bool_t HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_;
  Bool_t HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_;
  Bool_t HLT_Rsq0p35_;
  Bool_t HLT_Rsq0p40_;
  Bool_t HLT_RsqMR300_Rsq0p09_MR200_;
  Bool_t HLT_RsqMR320_Rsq0p09_MR200_;
  Bool_t HLT_RsqMR300_Rsq0p09_MR200_4jet_;
  Bool_t HLT_RsqMR320_Rsq0p09_MR200_4jet_;
  Bool_t HLT_IsoMu27_MET90_;
  Bool_t HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_;
  Bool_t HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_;
  Bool_t HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_;
  Bool_t HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_;
  Bool_t HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_;
  Bool_t HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_;
  Bool_t HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_;
  Bool_t HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_;
  Bool_t HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_;
  Bool_t HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_;
  Bool_t HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_;
  Bool_t HLT_PFMET100_PFMHT100_IDTight_PFHT60_;
  Bool_t HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_;
  Bool_t HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_;
  Bool_t HLT_Mu18_Mu9_SameSign_;
  Bool_t HLT_Mu18_Mu9_SameSign_DZ_;
  Bool_t HLT_Mu18_Mu9_;
  Bool_t HLT_Mu18_Mu9_DZ_;
  Bool_t HLT_Mu20_Mu10_SameSign_;
  Bool_t HLT_Mu20_Mu10_SameSign_DZ_;
  Bool_t HLT_Mu20_Mu10_;
  Bool_t HLT_Mu20_Mu10_DZ_;
  Bool_t HLT_Mu23_Mu12_SameSign_;
  Bool_t HLT_Mu23_Mu12_SameSign_DZ_;
  Bool_t HLT_Mu23_Mu12_;
  Bool_t HLT_Mu23_Mu12_DZ_;
  Bool_t HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_;
  Bool_t HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_;
  Bool_t HLT_DoubleMu3_DCA_PFMET50_PFMHT60_;
  Bool_t HLT_TripleMu_5_3_3_Mass3p8_DCA_;
  Bool_t HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;
  Bool_t HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;
  Bool_t HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;
  Bool_t HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_;
  Bool_t HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_;
  Bool_t HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_;
  Bool_t HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_;
  Bool_t HLT_QuadPFJet98_83_71_15_;
  Bool_t HLT_QuadPFJet103_88_75_15_;
  Bool_t HLT_QuadPFJet105_88_76_15_;
  Bool_t HLT_QuadPFJet111_90_80_15_;
  Bool_t HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_;
  Bool_t HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_;
  Bool_t HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_;
  Bool_t HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_;
  Bool_t HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_;
  Bool_t HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55_;
  Bool_t HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_;
  Bool_t HLT_Mu12_IP6_part0_;
  Bool_t HLT_Mu12_IP6_part1_;
  Bool_t HLT_Mu12_IP6_part2_;
  Bool_t HLT_Mu12_IP6_part3_;
  Bool_t HLT_Mu12_IP6_part4_;
  Bool_t HLT_Mu9_IP5_part0_;
  Bool_t HLT_Mu9_IP5_part1_;
  Bool_t HLT_Mu9_IP5_part2_;
  Bool_t HLT_Mu9_IP5_part3_;
  Bool_t HLT_Mu9_IP5_part4_;
  Bool_t HLT_Mu7_IP4_part0_;
  Bool_t HLT_Mu7_IP4_part1_;
  Bool_t HLT_Mu7_IP4_part2_;
  Bool_t HLT_Mu7_IP4_part3_;
  Bool_t HLT_Mu7_IP4_part4_;
  Bool_t HLT_Mu9_IP4_part0_;
  Bool_t HLT_Mu9_IP4_part1_;
  Bool_t HLT_Mu9_IP4_part2_;
  Bool_t HLT_Mu9_IP4_part3_;
  Bool_t HLT_Mu9_IP4_part4_;
  Bool_t HLT_Mu8_IP5_part0_;
  Bool_t HLT_Mu8_IP5_part1_;
  Bool_t HLT_Mu8_IP5_part2_;
  Bool_t HLT_Mu8_IP5_part3_;
  Bool_t HLT_Mu8_IP5_part4_;
  Bool_t HLT_Mu8_IP6_part0_;
  Bool_t HLT_Mu8_IP6_part1_;
  Bool_t HLT_Mu8_IP6_part2_;
  Bool_t HLT_Mu8_IP6_part3_;
  Bool_t HLT_Mu8_IP6_part4_;
  Bool_t HLT_Mu9_IP6_part0_;
  Bool_t HLT_Mu9_IP6_part1_;
  Bool_t HLT_Mu9_IP6_part2_;
  Bool_t HLT_Mu9_IP6_part3_;
  Bool_t HLT_Mu9_IP6_part4_;
  Bool_t HLT_Mu8_IP3_part0_;
  Bool_t HLT_Mu8_IP3_part1_;
  Bool_t HLT_Mu8_IP3_part2_;
  Bool_t HLT_Mu8_IP3_part3_;
  Bool_t HLT_Mu8_IP3_part4_;
  Bool_t HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_;
  Bool_t HLT_TrkMu6NoFiltersNoVtx_;
  Bool_t HLT_TrkMu16NoFiltersNoVtx_;
  Bool_t HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_;
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
  Bool_t Flag_ecalBadCalibFilter_;
  Bool_t Flag_goodVertices_;
  Bool_t Flag_eeBadScFilter_;
  Bool_t Flag_ecalLaserCorrFilter_;
  Bool_t Flag_trkPOGFilters_;
  Bool_t Flag_chargedHadronTrackResolutionFilter_;
  Bool_t Flag_muonBadTrackFilter_;
  Bool_t Flag_BadChargedCandidateFilter_;
  Bool_t Flag_BadPFMuonFilter_;
  Bool_t Flag_BadChargedCandidateSummer16Filter_;
  Bool_t Flag_BadPFMuonSummer16Filter_;
  Bool_t Flag_trkPOG_manystripclus53X_;
  Bool_t Flag_trkPOG_toomanystripclus53X_;
  Bool_t Flag_trkPOG_logErrorTooManyClusters_;
  Bool_t Flag_METFilters_;
  Bool_t L1Reco_step_;
  bool are_Jet_loaded_;
  vector<Jets> Jet_;
  bool are_IsoTrack_loaded_;
  vector<Isotracks> IsoTrack_;
  bool are_GenJetAK8_loaded_;
  vector<Genjetak8s> GenJetAK8_;
  bool are_GenVisTau_loaded_;
  vector<Genvistaus> GenVisTau_;
  bool are_CaloMET_loaded_;
  Calomet CaloMET_;
  bool are_GenDressedLepton_loaded_;
  vector<Gendressedleptons> GenDressedLepton_;
  bool are_LHEPart_loaded_;
  vector<Lheparts> LHEPart_;
  bool are_PV_loaded_;
  Pv PV_;
  bool are_Generator_loaded_;
  Generator Generator_;
  bool are_TrigObj_loaded_;
  vector<Trigobjs> TrigObj_;
  bool are_Photon_loaded_;
  vector<Photons> Photon_;
  bool are_GenJet_loaded_;
  vector<Genjets> GenJet_;
  bool are_RawMET_loaded_;
  Rawmet RawMET_;
  bool are_btagWeight_loaded_;
  Btagweight btagWeight_;
  bool are_SoftActivityJet_loaded_;
  Softactivityjet SoftActivityJet_;
  bool are_L1simulation_loaded_;
  L1Simulation L1simulation_;
  bool are_GenPart_loaded_;
  vector<Genparts> GenPart_;
  bool are_LHE_loaded_;
  Lhe LHE_;
  bool are_TkMET_loaded_;
  Tkmet TkMET_;
  bool are_Tau_loaded_;
  vector<Taus> Tau_;
  bool are_Flag_loaded_;
  Flag Flag_;
  bool are_PuppiMET_loaded_;
  Puppimet PuppiMET_;
  bool are_Muon_loaded_;
  vector<Muons> Muon_;
  bool are_L1Reco_loaded_;
  L1Reco L1Reco_;
  bool are_OtherPV_loaded_;
  vector<Otherpvs> OtherPV_;
  bool are_HLT_loaded_;
  Hlt HLT_;
  bool are_SV_loaded_;
  vector<Svs> SV_;
  bool are_Electron_loaded_;
  vector<Electrons> Electron_;
  bool are_MET_loaded_;
  Met MET_;
  bool are_GenMET_loaded_;
  Genmet GenMET_;
  bool are_FatJet_loaded_;
  vector<Fatjets> FatJet_;
  bool are_SubJet_loaded_;
  vector<Subjets> SubJet_;
  bool are_SubGenJetAK8_loaded_;
  vector<Subgenjetak8s> SubGenJetAK8_;
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

