#ifndef GENOBJECT_H
#define GENOBJECT_H

#include "URStreamer.h"
#include "PDGID.h"
#include <iostream>

class GenObject : public TLorentzVector
{
private:
	int pdgId_;
	int status_;

public:
	GenObject(const Genparticle& gp) : TLorentzVector(gp), pdgId_(gp.pdgId()), status_(gp.status())
	{
	}
	GenObject(const Pst& gp) : TLorentzVector(gp), pdgId_(gp.pdgId()), status_(gp.status())
	{
	}
	GenObject(const Genjet& jet) : TLorentzVector(jet), pdgId_(0), status_(0)
	{
	}
	GenObject() : TLorentzVector(), pdgId_(0), status_(0)
	{
	}
	
	int pdgId() const {return pdgId_;}
	int status() const {return status_;}

  friend std::ostream & operator<<(std::ostream &os, const GenObject &obj);
};

class GenW { //for some reason the gen-level W it is not tracked
public:
  enum DecayType {INVALID=0, LEPTONIC=1, HADRONIC=2};
  GenObject *first, *second;
  DecayType type = INVALID;

  GenW():
    first(0), 
    second(0),
    type(DecayType::INVALID) {}

  GenW(GenObject *first_, GenObject *second_=0) {
    if(ura::PDGID::e <= abs(first_->pdgId()) && abs(first_->pdgId()) <= ura::PDGID::nu_tau) {
      type = LEPTONIC;
      if(abs(first_->pdgId()) % 2 == 1) { //charged lepton first
        first = first_;
        second = second_;
      }
      else {
        first = second_;
        second = first_;
      }
    }
    else {
      type = HADRONIC;
      if(first_->Pt() > second_->Pt()){
        first =first_;
        second = second_;
      }
      else {
        first = second_;
        second = first_;
      }
    }
  }

  GenObject *up()   {return (first->pdgId() % 2 == 0) ? first : second;}
  GenObject *down() {return (first->pdgId() % 2 == 0) ? second : first;}

  friend std::ostream & operator<<(std::ostream &os, const GenW& w);
};

class GenTop: public GenObject {
public:
  GenObject *b;
  GenW W;

  GenTop():
    GenObject(),
    b(0),
    W() {}

  GenTop(const GenObject& t, GenObject *b_, GenW & W_):
    GenObject(t),
    b(b_),
    W(W_) {}

  friend std::ostream & operator<<(std::ostream &os, const GenTop& t);
};

class GenTTBar: public TLorentzVector {
public:
  enum DecayType {INVALID=0, FULLHAD=1, SEMILEP=2, FULLLEP=3};
  GenTop top, tbar;
  DecayType type;

  GenTTBar():
    TLorentzVector(),
    top(),
    tbar(),
    type(DecayType::INVALID) {}

  GenTTBar(const GenTop &t, const GenTop &topbar):
    TLorentzVector(t+topbar),
    top(t),
    tbar(topbar) {
    if(top.W.type == GenW::DecayType::INVALID || tbar.W.type == GenW::DecayType::INVALID)
      type = DecayType::INVALID;
    else if(top.W.type == GenW::DecayType::HADRONIC && tbar.W.type == GenW::DecayType::HADRONIC)
      type = DecayType::FULLHAD;
    else if(top.W.type == GenW::DecayType::LEPTONIC && tbar.W.type == GenW::DecayType::LEPTONIC)
      type = DecayType::FULLLEP;
    else 
      type = DecayType::SEMILEP;
  }

  static GenTTBar from_collections( vector<GenObject*>& wpartons, vector<GenObject*>& charged_leps,
            vector<GenObject*>& neutral_leps, GenObject* b, GenObject* bbar,
            GenObject* top, GenObject* tbar);

  GenTop* lep_top() {
    if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLHAD) return 0;
    return (top.W.type == GenW::DecayType::LEPTONIC) ? &top : &tbar;
  }

  GenTop* had_top() {
    if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLLEP) return 0;
    return (top.W.type == GenW::DecayType::HADRONIC) ? &top : &tbar;
  }

  GenObject* lepton(){
    if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLHAD) return 0;
    return (top.W.type == GenW::DecayType::LEPTONIC) ? top.W.first : tbar.W.first;
  }

  GenObject* lep_b() {
    if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLHAD) return 0;
    return (top.W.type == GenW::DecayType::LEPTONIC) ? top.b : tbar.b;
  }

  GenObject* had_b() {
    if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLLEP) return 0;
    return (top.W.type == GenW::DecayType::HADRONIC) ? top.b : tbar.b;
  }

  GenW* had_W() {
    if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLLEP) return 0;
    return (top.W.type == GenW::DecayType::HADRONIC) ? &top.W : &tbar.W;
  }

  friend std::ostream & operator<<(std::ostream &os, const GenTTBar& tt);
};

#endif
