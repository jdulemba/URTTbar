#include "GenObject.h"
#include "Logger.h"

std::ostream & operator<<(std::ostream &os, const GenObject &obj) {
  return os << "id: " << obj.pdgId_ << ", (pt, eta): (" << obj.Pt() << ", " << obj.Eta() << ")";
}

std::ostream & operator<<(std::ostream &os, const GenW& w) {
  os << "type: " << w.type;
  if(w.first) os << " --> " << *w.first << std::endl;
  else        os << " --> (NULL)" << std::endl;
  if(w.second) os << "        --> " <<  *w.second;
  else         os << "        --> (NULL)";
  return os;
}

std::ostream & operator<<(std::ostream &os, const GenTop& t) {
  if(t.b) os << "--> b: " << *t.b  << std::endl;
  else    os << "--> b: (NULL)" << std::endl;
  return os << "--> W: " << t.W;
}
std::ostream & operator<<(std::ostream &os, const GenTTBar& tt) {
  return os << "tt type: " << tt.type << std::endl
            << "top  " << tt.top << std::endl 
            << "tbar " << tt.tbar;
}

GenTTBar GenTTBar::from_collections( vector<GenObject*>& wpartons, vector<GenObject*>& charged_leps,
                    vector<GenObject*>& neutral_leps, GenObject* b, GenObject* bbar,
                    GenObject* top, GenObject* tbar) {
  GenW wplus, wminus;

  //leptonic Ws
  for(auto &lepton : charged_leps) {
    int lep_id = lepton->pdgId();
    int nu_id = -1*(lep_id + TMath::Abs(lep_id)/lep_id);
    for(auto &neutrino : neutral_leps) {
      if(neutrino->pdgId() == nu_id) {
        if(lep_id > 0) { //Wminus
          wminus = GenW(lepton, neutrino);
        } else {
          wplus = GenW(lepton, neutrino);            
        }
      }
    }
  } // for(auto &lepton : charged_leps_) {
  
    //hadronic Ws
  for(size_t p1_idx = 0; p1_idx < wpartons.size(); ++p1_idx) {
    for(size_t p2_idx = (p1_idx+1); p2_idx < wpartons.size(); ++p2_idx) {
      int p1_id = wpartons[p1_idx]->pdgId();
      int p2_id = wpartons[p2_idx]->pdgId();
      //check for u-d pairing
      if((TMath::Abs(p1_id) % 2) == (TMath::Abs(p2_id) % 2)) continue;
      //check for u-dbar or dbar-u
      if(p1_id*p2_id > 0) continue;
      //find up one
      int up_id = (TMath::Abs(p1_id) % 2 == 0)? p1_id : p2_id;
      if(up_id > 0) {
        wplus = GenW(wpartons[p1_idx], wpartons[p2_idx]);
      } else {
        wminus = GenW(wpartons[p1_idx], wpartons[p2_idx]);
      }
    }      
  } //for(size_t p1_idx = 0; p1_idx < wpartons_.size(); ++p1_idx)

  return GenTTBar(
    GenTop(*top, b, wplus), //top
    GenTop(*tbar, bbar, wminus)  //tbar
    );
}
