#ifndef Hypotheses_h
#define Hypotheses_h

#include "TLorentzVector.h"

class Permutation;
class GenTTBar;

namespace hyp {
  class Decay: public TLorentzVector {
  public:
    Decay():
      TLorentzVector() {}
    Decay(double x, double y, double z, double t):
      TLorentzVector(x, y, z, t) {}
    Decay(const TLorentzVector &lv):
      TLorentzVector(lv) {}
  
    void boost(const TVector3 &v);
    TVector3 unit3D() { return this->Vect().Unit();}
    void setv(TLorentzVector v) {this->SetPxPyPzE(v.Px(), v.Py(), v.Pz(), v.E());}
    double decay_opening_cm();
    TVector3 decay_plane();
    friend std::ostream & operator<<(std::ostream &os, const Decay &obj);
  private:
    virtual Decay* fst() = 0;
    virtual Decay* snd() = 0;
  };

  class NoDecay: public Decay {
  public:
    NoDecay():
      Decay() {}
    NoDecay(double x, double y, double z, double t):
      Decay(x, y, z, t) {}
    NoDecay(const TLorentzVector &lv):
      Decay(lv) {}
  private:
    virtual Decay* fst() {return NULL;}
    virtual Decay* snd() {return NULL;}
  };

  class WHyp: public Decay {
  private:
    NoDecay up_, dw_;
    int charge_=0;

    virtual Decay* fst() {return &up_;}
    virtual Decay* snd() {return &dw_;}

  public:
    WHyp():
      Decay(),
      up_(),
      dw_() {}

    WHyp to_CM();

    NoDecay& up() {return up_;}
    NoDecay& down() {return dw_;}
    NoDecay& l() {return dw_;}
    NoDecay& nu() {return up_;}
    void up(TLorentzVector up) {up_ = up;}
    void down(TLorentzVector dw) {dw_ = dw;}
    void l(TLorentzVector v) {dw_ = v;}
    void nu(TLorentzVector v) {up_ = v;}
  };

  class Top: public Decay {
  private:
    NoDecay b_;
    WHyp w_;

    virtual Decay* fst() {return &b_;}
    virtual Decay* snd() {return &w_;}

  public:
    Top to_CM();

    //constructors
    Top():
      Decay(),
      b_(),
      w_() {}

    //element access
    NoDecay& b() {return b_;}
    WHyp& W() {return w_;}
    void b(TLorentzVector b) {b_ = b;}
    void W(WHyp w) {w_ = w;}
  };

  class TTbar: public Decay {
  private:
    Top t_, tbar_;
    bool t_leptonic_=false;

    virtual Decay* fst() {return &t_;}
    virtual Decay* snd() {return &tbar_;}

  public:
    TTbar to_CM();
  
    //constructors
    TTbar(Permutation &p);
    TTbar(GenTTBar &g);
    TTbar():
      Decay(),
      t_leptonic_(false),
      t_(),
      tbar_() {}

    //element access
    Top& top()  {return t_;}
    Top& tbar() {return tbar_;}
    Top& tlep() {return (t_leptonic_) ? t_    : tbar_;}
    Top& thad() {return (t_leptonic_) ? tbar_ : t_   ;}
    int lep_charge() {return (t_leptonic_) ? 1 : -1;}

    void top( Top t) {t_ = t;}
    void tbar(Top t) {tbar_ = t;}
  };
}//namespace


#endif
