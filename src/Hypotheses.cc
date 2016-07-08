#include "Analyses/URTTbar/interface/Hypotheses.h"
#include "Analyses/URTTbar/interface/TwoBodyDecay.h"
#include "Analyses/URTTbar/interface/Permutation.h"
#include "Analyses/URTTbar/interface/GenObject.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include <iostream>

namespace hyp {
//
// Decay
//
  typedef TLorentzVector tlv;
  void Decay::boost(const TVector3 &v) {
    this->Boost(v);
    if(fst()) fst()->boost(v);
    if(snd()) snd()->boost(v);
  }

  double Decay::decay_opening_cm() {
    if(!fst() || !snd()) return -2;
    TVector3 bv = this->BoostVector()*-1;
    TLorentzVector f(*fst()); f.Boost(bv);
    TLorentzVector s(*snd()); s.Boost(bv);
    return f.Dot(s);
  }
  TVector3 Decay::decay_plane() {
    if(!fst() || !snd()) return TVector3();
    return (fst()->unit3D().Cross(snd()->unit3D())).Unit();
  }

  std::ostream & operator<<(std::ostream &os, const Decay &obj) {
    return os << "Decay{" << obj.Px() << ", " << obj.Py() << ", " << obj.Pz() << ", " << obj.E() << "}";
  }
//
// TTbar
//
  TTbar TTbar::to_CM() {
    TTbar ret=*this; //this time no ptrs, it makes a hard clone
    ret.boost(-1*this->BoostVector());
    return ret;
  }

  TTbar::TTbar(Permutation &p):
    TTbar()
  {
    int lc = p.LepCharge();
    if(lc == 0) {
      Logger::log().fatal() << "Permutation lepton charge not properly set!" << std::endl;
      throw 42;
    }
    t_leptonic_ = (lc > 0); //lepton charge == +1 --> top leptonic
    tlep().b(*p.BLep());
    tlep().W().l(*p.L());
    tlep().W().nu(p.Nu());
    tlep().W().setv(*p.L()+p.Nu());
    tlep().setv(tlep().W()+*p.BLep());

    thad().b(*p.BHad());
		//if gen matching exists and makes sense
		if(p.WJa()->match() && p.WJb()->match() &&
			 (p.WJa()->match()->pdgId()+p.WJb()->match()->pdgId()) % 2 == 1) { //check that one is up and one is down
			if(p.WJa()->match()->pdgId() % 2 == 0) {
				thad().W().up(*p.WJa());
				thad().W().down(*p.WJb());				
			}
			else {
				thad().W().up(*p.WJb());
				thad().W().down(*p.WJa());				
			}
		} 
		else {
			thad().W().up(*p.WJa()); //Assumes WJa is up one!
			thad().W().down(*p.WJb());
		}
    thad().W().setv(*p.WJa()+*p.WJb());
    thad().setv(thad().W()+*p.BHad());
    setv(thad()+tlep());
  }

  TTbar::TTbar(GenTTBar &g):
    TTbar()
  {
    if(g.type == GenTTBar::DecayType::INVALID) {
      Logger::log().fatal() << "GenTTBar status is invalid!" << std::endl;
      throw 42;
    }

    t_leptonic_ = (g.top.W.type == GenW::DecayType::LEPTONIC);

    top().b(*g.top.b);
    top().W().up(*g.top.W.up());
    top().W().down(*g.top.W.down());
    top().W().setv(*g.top.W.up()+*g.top.W.down());
    top().setv(top().W()+top().b());

    tbar().b(*g.tbar.b);
    tbar().W().up(*g.tbar.W.up());
    tbar().W().down(*g.tbar.W.down());
    tbar().W().setv(*g.tbar.W.up()+*g.tbar.W.down());
    tbar().setv(tbar().W()+tbar().b());
    setv(top()+tbar());
  }

//
// Top
//
  Top Top::to_CM() {
    Top ret=*this; //this time no ptrs, it makes a hard clone
    ret.boost(this->BoostVector()*-1);
    return ret;
  }

//
// W
//
  WHyp WHyp::to_CM() {
    WHyp ret=*this; //this time no ptrs, it makes a hard clone
    ret.boost(this->BoostVector()*-1);
    return ret;
  }
}//namespace hyp 
