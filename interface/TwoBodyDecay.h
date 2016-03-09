#ifndef TwoBodyDecay_h
#define TwoBodyDecay_h

#include "TLorentzVector.h"
#include <memory>
#include <vector>
/*
Tree-like structure of a two-body decay chain
 */

class Permutation;
class GenTTBar;

class TwoBodyDecay: public TLorentzVector{
public:
  TwoBodyDecay(Permutation &p);
  TwoBodyDecay(GenTTBar &g);

  TwoBodyDecay(): //default
    TLorentzVector(),
    fst_(),
    snd_() {}

  TwoBodyDecay(TwoBodyDecay *one, TwoBodyDecay *two): //from decay product
    TLorentzVector(*one+*two),
    fst_(one),
    snd_(two) {}

  TwoBodyDecay(const TLorentzVector *one, const TLorentzVector *two): //from decay product
    TLorentzVector(*one+*two),
    fst_(new TwoBodyDecay(one)),
    snd_(new TwoBodyDecay(two)) {}

  TwoBodyDecay(const TLorentzVector *val): //last item
    TLorentzVector(*val),
    fst_(),
    snd_() {}

  TwoBodyDecay(const TLorentzVector *val, TwoBodyDecay *one, TwoBodyDecay *two): //decay products and total
    TLorentzVector(*val),
    fst_(one),
    snd_(two) {}

  TwoBodyDecay(const TLorentzVector *val, std::shared_ptr<TwoBodyDecay> one, std::shared_ptr<TwoBodyDecay> two): //decay products and total
    TLorentzVector(*val),
    fst_(one),
    snd_(two) {}

  std::shared_ptr<TwoBodyDecay> fst() {return fst_;}
  std::shared_ptr<TwoBodyDecay> snd() {return snd_;}

  std::shared_ptr<TwoBodyDecay> clone();

  void boost(const TVector3 &v);
  std::shared_ptr<TwoBodyDecay> to_CM();

  TVector3 unit3D() { return this->Vect().Unit(); }
  double cosThetaStar_decay();

private:
  std::shared_ptr<TwoBodyDecay> fst_;
  std::shared_ptr<TwoBodyDecay> snd_;
};

#endif
