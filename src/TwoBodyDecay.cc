#include "TwoBodyDecay.h"
#include "Permutation.h"
#include "GenObject.h"

using namespace std;

TwoBodyDecay::TwoBodyDecay(Permutation &p):
  TwoBodyDecay(
    new TwoBodyDecay(
      new TwoBodyDecay(p.BLep()),
      new TwoBodyDecay(p.L(), p.NuPtr())
      ),
    new TwoBodyDecay(
      new TwoBodyDecay(p.BHad()),
      new TwoBodyDecay(p.WJa(), p.WJb())
      )
    ) {}

TwoBodyDecay::TwoBodyDecay(GenTTBar &g):
  TwoBodyDecay(
    new TwoBodyDecay(
      &g.top,
      new TwoBodyDecay(g.top.b),
      new TwoBodyDecay(g.top.W.first, g.top.W.second)
      ),
    new TwoBodyDecay(
      &g.tbar,
      new TwoBodyDecay(g.tbar.b),
      new TwoBodyDecay(g.tbar.W.first, g.top.W.second)
      )
    )  {}       

std::shared_ptr<TwoBodyDecay> TwoBodyDecay::clone() {
  shared_ptr<TwoBodyDecay> fst = (fst_.get()) ? fst_->clone() : 0;
  shared_ptr<TwoBodyDecay> snd = (snd_.get()) ? snd_->clone() : 0;  
  std::shared_ptr<TwoBodyDecay> ret(new TwoBodyDecay(this, fst, snd));
  return ret;
}

void TwoBodyDecay::boost(const TVector3 &v) {
  this->Boost(v);
  if(fst_.get()) fst_->Boost(v);
  if(snd_.get()) snd_->Boost(v);
}

std::shared_ptr<TwoBodyDecay> TwoBodyDecay::to_CM() {
  shared_ptr<TwoBodyDecay> ret = this->clone();
  ret->boost(this->BoostVector()*-1);
  return ret;
}

double TwoBodyDecay::cosThetaStar_decay() {
  if(!fst().get()) return -2;
  TVector3 direction = unit3D();
  shared_ptr<TwoBodyDecay> cmframe = this->to_CM();
  TVector3 fst_dir = cmframe->fst()->unit3D();
  return direction.Dot(fst_dir);
}
// 
// 
// 
// 
// 
// 
// 
// 
// 
