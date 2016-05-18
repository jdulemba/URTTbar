#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include "TMath.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "URAnalysis/AnalysisFW/interface/PDGID.h"

using namespace TMath;

std::ostream & operator<<(std::ostream &os, const Genparticle& w) {
  os << "Genparticle(" << w.idx() << ", id:";
  if(ura::pdg_names.find(w.pdgId()) != ura::pdg_names.end()){
    os << ura::pdg_names.at(w.pdgId());// << ", stat:" << w.status() << ")";
  } else {
    os << w.pdgId();
  }    
  return os << ", stat:" << w.status() << ")";//, " 
  //<< w.Px() << ", " << w.Py() << ", " << w.Pz() << ", " << w.E() << ")";
}

TTGenParticleSelector::TTGenParticleSelector(SelMode mode):
  lhes_(),
  selected_(),
  wpartons_(),
  charged_leps_(),
  neutral_leps_(),
  final_charged_leps_(),
  added_jets_(),
  jets_(),  
  top_(0),
  tbar_(0),
  b_(0),
  bbar_(0),
  topcounter_(0),
  lepdecays_(0),
  ttbar_(),
  ttbar_final_() {
  setmode(mode);

  URParser &parser = URParser::instance();
  parser.addCfgParameter<float>("gen_jets", "ptmin", "minimum pt");
  parser.addCfgParameter<float>("gen_jets", "etamax", "maximum eta");

  parser.addCfgParameter<float>("acceptance", "lepton_ptmin", "minimum lepton pt", 0.);
  parser.addCfgParameter<float>("acceptance", "lepton_etamax", "maximum lepton eta", 999.);

  parser.addCfgParameter<float>("acceptance", "wjet_ptsoft", "minimum W-jet pt", 0.);
  parser.addCfgParameter<float>("acceptance", "wjet_pthard", "minimum leading W-jet pt", 0.);
  parser.addCfgParameter<float>("acceptance", "wjet_etamax", "maximum W-jet eta", 999.);

  parser.addCfgParameter<float>("acceptance", "bjet_ptsoft", "minimum b-jet pt", 0.);
  parser.addCfgParameter<float>("acceptance", "bjet_pthard", "minimum leading b-jet pt", 0.);
  parser.addCfgParameter<float>("acceptance", "bjet_etamax", "maximum b-jet eta", 999.);

  parser.addCfgParameter<float>("acceptance", "parton_separation", "parton separation");

  parser.parseArguments();
  parser.parseArguments();

  jetptmin_  = parser.getCfgPar<float>("gen_jets", "ptmin" );
  jetetamax_ = parser.getCfgPar<float>("gen_jets", "etamax");

  cut_lep_ptmin_ = parser.getCfgPar<float>("acceptance", "lepton_ptmin"  );
  cut_lep_etamax_ = parser.getCfgPar<float>("acceptance", "lepton_etamax");
  cut_wjet_ptsoft_ = parser.getCfgPar<float>("acceptance", "wjet_ptsoft" );
  cut_wjet_pthard_ = parser.getCfgPar<float>("acceptance", "wjet_pthard" );
  cut_wjet_etamax_ = parser.getCfgPar<float>("acceptance", "wjet_etamax" );
  cut_bjet_ptsoft_ = parser.getCfgPar<float>("acceptance", "bjet_ptsoft" );
  cut_bjet_pthard_ = parser.getCfgPar<float>("acceptance", "bjet_pthard" );
  cut_bjet_etamax_ = parser.getCfgPar<float>("acceptance", "bjet_etamax" );
  cut_parton_separation_ = parser.getCfgPar<float>("acceptance", "parton_separation");
}
 

void TTGenParticleSelector::select_pstop(URStreamer& event) 
{
  const vector<Pst>& pseudotops = event.PSTs();
  if(pseudotops.size() == 10) {
    topcounter_ = 2;

    if(Abs(pseudotops[8].pdgId()) == 11 || Abs(pseudotops[8].pdgId()) == 13) {
      selected_.push_back(pseudotops[8]);
      charged_leps_.push_back(&(selected_.back()));
      final_charged_leps_.push_back(&(selected_.back()));
      selected_.push_back(pseudotops[9]);
      neutral_leps_.push_back(&(selected_.back()));
      lepdecays_++;
    }
    else {
      selected_.push_back(pseudotops[8]);
      wpartons_.push_back(&(selected_.back()));
      selected_.push_back(pseudotops[9]);
      wpartons_.push_back(&(selected_.back()));
    }

    if(Abs(pseudotops[3].pdgId()) == 11 || Abs(pseudotops[3].pdgId()) == 13) {
      selected_.push_back(pseudotops[3]);
      charged_leps_.push_back(&(selected_.back()));
      final_charged_leps_.push_back(&(selected_.back()));
      selected_.push_back(pseudotops[4]);
      neutral_leps_.push_back(&(selected_.back()));
      lepdecays_++;
    }
    else {
      selected_.push_back(pseudotops[3]);
      wpartons_.push_back(&(selected_.back()));
      selected_.push_back(pseudotops[4]);
      wpartons_.push_back(&(selected_.back()));
    }

    if(pseudotops[2].pdgId() == 5) {
      selected_.push_back(pseudotops[2]);
      b_ = &(selected_.back());
      selected_.push_back(pseudotops[7]);
      bbar_ = &(selected_.back());
    }
    else {
      selected_.push_back(pseudotops[2]);
      bbar_ = &(selected_.back());
      selected_.push_back(pseudotops[7]);
      b_ = &(selected_.back());
    }
  }
}

void TTGenParticleSelector::select_herwig(URStreamer& event)
{
  const vector<Genparticle>& gps = event.genParticles();
  for(vector<Genparticle>::const_iterator gp = gps.begin(); gp != gps.end(); ++gp) {
    if(gp->status() == 11) {
      if(gp->pdgId() == 6 && gps[gp->firstDaughtIdx()].pdgId() != 6) {
        //weight *= 1.+ctopptweight*(gp->Pt()-200.)/1000.;
        topcounter_++;
        selected_.push_back(*gp);
        top_ = &(selected_.back());
      }
      else if(gp->pdgId() == -6 && gps[gp->firstDaughtIdx()].pdgId() != -6) {
        topcounter_++;
        selected_.push_back(*gp);
        tbar_ = &(selected_.back());
      }
      else if(gp->pdgId() == 5 && gps[gp->momIdx()[0]].pdgId() == 6) {
        selected_.push_back(*gp);
        b_ = &(selected_.back());
      }
      else if(gp->pdgId() == -5 && gps[gp->momIdx()[0]].pdgId() == -6) {
        selected_.push_back(*gp);
        bbar_ = &(selected_.back());
      }
      else if(Abs(gp->pdgId()) < 6 && gp->momIdx().size() != 0 && Abs(gps[gp->momIdx()[0]].pdgId()) == 24) {
        selected_.push_back(*gp);
        wpartons_.push_back(&(selected_.back()));
      }
      else if(gp->momIdx().size() != 0 && Abs(gps[gp->momIdx()[0]].pdgId()) == 24) {
        if(Abs(gp->pdgId()) == 11 || Abs(gp->pdgId()) == 13) {
          selected_.push_back(*gp);
          charged_leps_.push_back(&(selected_.back()));
          final_charged_leps_.push_back(&(selected_.back()));	
        }
        if(Abs(gp->pdgId()) == 12 || Abs(gp->pdgId()) == 14) {
          selected_.push_back(*gp);
          neutral_leps_.push_back(&(selected_.back()));	
          lepdecays_++;
        }
        if(Abs(gp->pdgId()) == 16) {
          lepdecays_++;
        }
      }
    }
  }
}

void TTGenParticleSelector::select_normal(URStreamer& event)
{
  const vector<Genparticle>& gps = event.genParticles();
  for(vector<Genparticle>::const_iterator gp = gps.begin(); gp != gps.end(); ++gp) {
    if(gp->status() > 21 && gp->status() < 30 && gp->momIdx().size() != 0) {
      if(gp->pdgId() == 6) {
        //if(gps[gp->momIdx()[0]].pdgId() != 21) Logger::log().info() << gp->pdgId() << " <-- " << gps[gp->momIdx()[0]].pdgId() << std::endl;
        topcounter_++;
        selected_.push_back(*gp);
        top_ = &(selected_.back());
      }
      else if(gp->pdgId() == -6) {
        //if(gps[gp->momIdx()[0]].pdgId() != 21) Logger::log().info() << gp->pdgId() << " <-- " << gps[gp->momIdx()[0]].pdgId() << std::endl;
        topcounter_++;
        selected_.push_back(*gp);
        tbar_ = &(selected_.back());
      }
      else if(gp->pdgId() == 5 && gps[gp->momIdx()[0]].pdgId() != 24) {
        selected_.push_back(*gp);
        b_ = &(selected_.back());
      }
      else if(gp->pdgId() == -5 && gps[gp->momIdx()[0]].pdgId() != -24) {
        selected_.push_back(*gp);
        bbar_ = &(selected_.back());
      }
      else if(Abs(gp->pdgId()) < 6 && Abs(gps[gp->momIdx()[0]].pdgId()) == w_decay_momid_) {
        selected_.push_back(*gp);
        wpartons_.push_back(&(selected_.back()));
      }
    }

    if(gp->momIdx().size() != 0 && Abs(gps[gp->momIdx()[0]].pdgId()) == w_decay_momid_) {	
      if(Abs(gp->pdgId()) == 11 || Abs(gp->pdgId()) == 13) {
        selected_.push_back(*gp);
        charged_leps_.push_back(&(selected_.back()));
      }
      if(Abs(gp->pdgId()) == 12 || Abs(gp->pdgId()) == 14) {
        selected_.push_back(*gp);
        neutral_leps_.push_back(&(selected_.back()));	
        lepdecays_++;
      }
      if(Abs(gp->pdgId()) == 16) {
        lepdecays_++;
      }
    }

    if(gp->status() == 1 && gp->momIdx().size() != 0 && (Abs(gps[gp->momIdx()[0]].pdgId()) == 24 || gp->pdgId() == gps[gp->momIdx()[0]].pdgId())) {
      if(Abs(gp->pdgId()) == 11 || Abs(gp->pdgId()) == 13) {
        selected_.push_back(*gp);
        final_charged_leps_.push_back(&(selected_.back()));	
      }
    }
  }
}

int TTGenParticleSelector::comes_from_top(LHEParticle &lhe) {
  pair<int,int> moms = lhe.mothers_range();
  if(moms.first != moms.second || moms.first < 0) return -1; //not coming from top
  LHEParticle &mom = lhes_[moms.first];
  if(Abs(mom.pdgId()) == ura::PDGID::t) return moms.first;
  return comes_from_top(mom);
}

void TTGenParticleSelector::select_lhe(URStreamer& event) {
  lhes_ = LHEParticle::LHEParticles(event);
  for(auto &lhe : lhes_) {
    if(lhe.pdgId() == ura::PDGID::t) {
      selected_.push_back(lhe);
      topcounter_++;
      top_ = &(selected_.back());
      continue;
    }
    else if(lhe.pdgId() == ura::PDGID::tbar) {
      topcounter_++;
      selected_.push_back(lhe);
      tbar_ = &(selected_.back());
      continue;
    }
    if(lhe.status() != 1) {continue;}
    int top_id = this->comes_from_top(lhe);
    if(top_id == 0) {continue;}//does not come from top
    selected_.push_back(lhe);
    int mom = lhe.mothers_range().first;

    if(lhe.pdgId() == ura::PDGID::b && lhes_[mom].pdgId() != ura::PDGID::Wplus) {
      b_ = &(selected_.back());
    }
    else if(lhe.pdgId() == ura::PDGID::bbar && lhes_[mom].pdgId() != ura::PDGID::Wminus) {
      bbar_ = &(selected_.back());
    }
    else if(Abs(lhes_[mom].pdgId()) == w_decay_momid_) {
      if(Abs(lhe.pdgId()) < 6) wpartons_.push_back(&(selected_.back()));
      if(Abs(lhe.pdgId()) == ura::PDGID::e || Abs(lhe.pdgId()) == ura::PDGID::mu) charged_leps_.push_back(&(selected_.back()));
      if(Abs(lhe.pdgId()) == ura::PDGID::nu_e || Abs(lhe.pdgId()) == ura::PDGID::nu_mu) {
        neutral_leps_.push_back(&(selected_.back()));
        lepdecays_++;
      }
      if(Abs(lhe.pdgId()) == ura::PDGID::tau || Abs(lhe.pdgId()) == ura::PDGID::nu_tau) {
        lepdecays_++;
      }
    }
  }
}

void TTGenParticleSelector::reset() {
  //RESETS NEEDED VECTORS
  lhes_.clear();
  selected_.clear();
  wpartons_.clear();
  charged_leps_.clear();
  neutral_leps_.clear();
  final_charged_leps_.clear();
  top_ = 0;
  tbar_ = 0;
  b_ = 0;
  bbar_ = 0;
  // wplus_ = 0;
  // wminus_ = 0;
  topcounter_ = 0;
  ttbar_.type = GenTTBar::DecayType::INVALID;
  ttbar_final_.type = GenTTBar::DecayType::INVALID;
  jets_.clear();
  added_jets_.clear();
  is_in_acceptance_ =-1;
}

bool  TTGenParticleSelector::select(URStreamer& event) {
  reset();
  //SELECT BASED ON SELECTION MODE
  //TODO
  switch(mode_) {
  case NORMAL: 
  case MADGRAPH:  select_normal(event); break;
  case PSEUDOTOP: select_pstop(event); break;
  case HERWIGPP:  select_herwig(event); break;
  case LHE: 
  case MADGRAPHLHE: select_lhe(event); break;
  case FULLDEP:     select_with_deps(event); break;
  }

  //Build GenTTBar
  if(mode_ == PSEUDOTOP) {
    Logger::log().debug() << "pstops: " << event.PSTs().size() << " pst leps: " << event.PSTleptons().size() << " pst jets: " << event.PSTjets().size() << " pst nus: " << event.PSTneutrinos().size() << std::endl;
    Logger::log().debug() << "wpartons: " << wpartons_.size() << ", charged_leps: " << charged_leps_.size()
                          << ", neutral_leps: " <<neutral_leps_.size() << ", b: " << b_ << ", bbar_: "
                          << bbar_ << ", top: " << top_ << ", tbar: " << tbar_ << std::endl;
    if(wpartons_.size()+charged_leps_.size() == 0 || !b_ || !bbar_) {
      Logger::log().error() << "Wrong matching, returning" << std::endl;
      return false;
    }
  }
  ttbar_ = GenTTBar::from_collections( 
    wpartons_, charged_leps_, neutral_leps_, 
    b_, bbar_, top_, tbar_);
  ttbar_final_ = GenTTBar::from_collections(
    wpartons_, final_charged_leps_, neutral_leps_,
    b_, bbar_, top_, tbar_);
  // Logger::log().debug() << "ttbar_: " << ttbar_ <<std::endl;
  // Logger::log().debug() << "ttbar_final_: " << ttbar_final_ <<std::endl;

  if(!ttbar_.is_complete()) return false;

  //Makes collection of gen jets not in the partons
	if(ttbar_.type == GenTTBar::DecayType::SEMILEP) {
    if(!charged_leps_[0] || !wpartons_[0] || !wpartons_[1] || !bbar_ || !b_) return false;

		const vector<Genjet>& genjets = event.genjets();
		for(vector<Genjet>::const_iterator gj = genjets.begin(); gj != genjets.end(); ++gj)
		{
			if(gj->Pt() < jetptmin_ || Abs(gj->Eta()) > jetetamax_) {continue;}
			if(gj->DeltaR(*charged_leps_[0]) < 0.4) continue;
			if(gj->DeltaR(*wpartons_[0]) < 0.4)	continue;
			if(gj->DeltaR(*wpartons_[1]) < 0.4)	continue;
			if(gj->DeltaR(*bbar_) < 0.4) continue;
			if(gj->DeltaR(*b_) < 0.4) continue;

			jets_.push_back(Genjet(*gj));
			added_jets_.push_back(&(jets_.back()));
		}
	}
  return true;
}

bool TTGenParticleSelector::is_in_acceptance(GenTTBar::DecayType decay_mode) {
  if(is_in_acceptance_ > -1) return is_in_acceptance_; //caching!
  if(ttbar_.type != decay_mode) {
    is_in_acceptance_ = 0;
    return false;
  }
  if(decay_mode != GenTTBar::DecayType::SEMILEP) { //TODO
    Logger::log().error() << "The code is not yet ready to handle acceptances with decay modes different from semileptonic! Fix it!" << std::endl;
    return false;
  }
  
  //lepton check
  if(ttbar_.lepton()->Pt() < cut_lep_ptmin_ || Abs(ttbar_.lepton()->Eta()) > cut_lep_etamax_) {
    is_in_acceptance_ = 0;
    return false;
  }

  //hadronic W check
  if(ttbar_.had_W()->first->Pt() < cut_wjet_pthard_ || ttbar_.had_W()->second->Pt()< cut_wjet_ptsoft_ ||
     Abs(ttbar_.had_W()->first->Eta()) > cut_wjet_etamax_ || Abs(ttbar_.had_W()->second->Eta()) > cut_wjet_etamax_) {
    is_in_acceptance_ = 0;
    return false;
  }

  //bjets check   
  if(ttbar_.top.b->Pt() < cut_bjet_ptsoft_ || ttbar_.tbar.b->Pt() < cut_bjet_ptsoft_){
    is_in_acceptance_ = 0;
    return false;
  }
  if(ttbar_.top.b->Pt() < cut_bjet_pthard_ && ttbar_.tbar.b->Pt() < cut_bjet_pthard_){
    is_in_acceptance_ = 0;
    return false;
  }
  if(Abs(ttbar_.top.b->Eta()) < cut_bjet_etamax_ || Abs(ttbar_.tbar.b->Eta()) < cut_bjet_etamax_){
    is_in_acceptance_ = 0;
    return false;
  }

  //parton separation
  if(ttbar_.had_W()->first->DeltaR(*ttbar_.had_W()->second) < cut_parton_separation_ || 
     ttbar_.had_W()->first->DeltaR(*ttbar_.top.b) < cut_parton_separation_ ||
     ttbar_.had_W()->first->DeltaR(*ttbar_.tbar.b) < cut_parton_separation_ ||
     ttbar_.had_W()->second->DeltaR(*ttbar_.top.b) < cut_parton_separation_ ||
     ttbar_.had_W()->second->DeltaR(*ttbar_.tbar.b) < cut_parton_separation_ ||
     ttbar_.top.b->DeltaR(*ttbar_.tbar.b) < cut_parton_separation_ ) {
    is_in_acceptance_ = 0;
    return false;
  }

  is_in_acceptance_ = 1;
  return true;
}

//
// OLD CODE, might come useful as matching is different
//

bool TTGenParticleSelector::descends(const Genparticle* mother, const Genparticle* child) {
  for(int mom_idx : child->momIdx()) 
    if(mom_idx == mother->idx()) 
      return true;
  return false;
}

vector< pair<const Genparticle*, const Genparticle*> > TTGenParticleSelector::Collapse(
  vector<const Genparticle*> &roots, 
  vector<const Genparticle*> &particles) {
  //particles MUST NOT contain roots
  vector< pair<const Genparticle*, const Genparticle*> > retval;
  vector<bool> used(particles.size(), false);
  
  for(auto& root : roots) {
    const Genparticle* current = root;
    bool go_on=true;
    while(go_on){
      go_on=false;
      //loop over particles
      for(size_t pos=0; pos<particles.size(); pos++){
        //check if already used and has mothers
        if(particles[pos]->momIdx().size() == 0 || used[pos]) continue;

        //check if is a descendent
        bool descends_from = descends(current, particles[pos]);
        
        //if so use it
        if(descends_from) {
          used[pos] = true;
          go_on=true;
          current = particles[pos];
          break;
        }
      }
    }
    
    //make pair
    retval.push_back( make_pair(root, current) );
  }

	return retval;
}

bool TTGenParticleSelector::assign(const Genparticle& gp, const std::vector<Genparticle>& gps, 
            vector<const Genparticle*> &collection, vector<const Genparticle*> &roots, 
            ura::PDGID to_match) {
  if(gp.pdgId() == to_match){
    if(gp.momIdx().size() == 0 || gps[gp.momIdx()[0]].pdgId() != to_match) {
      roots.push_back(&gp);
    } else {
      collection.push_back(&gp);
    }
    return true;
  }
  else return false;
}

void TTGenParticleSelector::select_with_deps(URStreamer& event)
{
	const std::vector<Genparticle>& gps = event.genParticles();

	//find top and tbar (multiple due to radiation
	vector<const Genparticle*> tops;
	vector<const Genparticle*> tbars;  
  vector<const Genparticle*> Wpluses;
  vector<const Genparticle*> Wminuses;
  vector<const Genparticle*> bs;
  vector<const Genparticle*> bbars;

	vector<const Genparticle*> root_tops;
	vector<const Genparticle*> root_tbars;  
  vector<const Genparticle*> root_Wpluses;
  vector<const Genparticle*> root_Wminuses;

  //find tops, bs, and Ws
	for(auto& gp : gps) {
    if( assign(gp, gps, tops, root_tops, ura::PDGID::t) ) continue;
    else if(assign(gp, gps, tbars, root_tbars, ura::PDGID::tbar)) continue;
    else if(assign(gp, gps, Wpluses, root_Wpluses, ura::PDGID::Wplus)) continue;
    else if(assign(gp, gps, Wminuses, root_Wminuses, ura::PDGID::Wminus)) continue;
    else if(gp.pdgId() == ura::PDGID::b) bs.push_back(&gp);
    else if(gp.pdgId() == ura::PDGID::bbar) bbars.push_back(&gp);
	}

  //collapse them (compress the same paricle)
	auto collapsed_tops = Collapse(root_tops, tops);
	auto collapsed_tbars = Collapse(root_tbars, tbars);
	auto collapsed_Wpluses = Collapse(root_Wpluses, Wpluses);
	auto collapsed_Wminuses = Collapse(root_Wminuses, Wminuses);

  if(collapsed_tops.size() != 1 || collapsed_tbars.size() != 1) {
    Logger::log().error() << "Could not find the proper number of tops!" << endl;
    throw 42;
  }

  //store tops
  topcounter_ += 2;
  selected_.push_back(*(collapsed_tops[0].second));
  top_ = &(selected_.back());
  selected_.push_back(*(collapsed_tbars[0].second));
  tbar_ = &(selected_.back());
  
  //select and store b quarks
  for(auto &bcan : bs) {
    if(descends(collapsed_tops[0].second, bcan) ) {
      selected_.push_back(*bcan);
      b_ =  &(selected_.back());
      break;
    }
  }
  for(auto &bcan : bbars) {
    if(descends(collapsed_tbars[0].second, bcan) ) {
      selected_.push_back(*bcan);
      bbar_ =  &(selected_.back());
      break;
    }
  }
  
  //select Ws
  const Genparticle* wplus = 0;
  const Genparticle* wminus = 0;
  for(auto &ws : collapsed_Wpluses) {
    if(descends(collapsed_tops[0].second, ws.first)) {
      wplus = ws.second;
    }
  }
  for(auto &ws : collapsed_Wminuses) {
    if(descends(collapsed_tbars[0].second, ws.first)) {
      wminus = ws.second;
    }
  }

  if(!wplus || !wminus) {
    Logger::log().error() << "Could not find the Ws from top decay!" << endl;
    throw 42;
  }
  
  //cout << wplus << " " << wminus << endl;

	//look for W decay products
  vector<const Genparticle*> root_leps;
  for(auto& gp : gps) {
    int abs_pdgid = fabs(gp.pdgId());
    // if(abs_pdgid > 5 && abs_pdgid % 2 == 0) {
    //   cout << "neutrino found: mother " << gps[gp.momIdx()[0]].pdgId() << endl;
    // }

    if(descends(wplus, &gp) || descends(wminus, &gp)) {
      selected_.push_back(gp);
      if(abs_pdgid < 5) wpartons_.push_back(&(selected_.back()));//quark      
      else if(abs_pdgid % 2 == 0) neutral_leps_.push_back(&(selected_.back())); //neutrino
      else {//charged lepton
        root_leps.push_back(&gp);
        charged_leps_.push_back(&(selected_.back()));
      }
    } 
  }

  //look for lepton radiation
  if(root_leps.size() > 0) {
    vector<const Genparticle*> leps;
    for(auto& gp : gps) {
      for(auto& root : root_leps) {
        if(root->pdgId() == gp.pdgId() && root->idx() != gp.idx()) leps.push_back(&gp);
      }
    }//for(auto& gp : gps)
    
    auto collapsed_leps = Collapse(root_leps, leps);
    for(auto& lep : collapsed_leps) {
      selected_.push_back(*(lep.second));
      final_charged_leps_.push_back(&(selected_.back()));
    }
  } //if(root_leps.size() > 0)

}


