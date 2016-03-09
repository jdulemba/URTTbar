#include "TTGenParticleSelector.h"
#include "TMath.h"
#include "URParser.h"
#include "Logger.h"
#include "PDGID.h"

using namespace TMath;

TTGenParticleSelector::TTGenParticleSelector(SelMode mode):
  mode_(mode),
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
        topcounter_++;
        selected_.push_back(*gp);
        top_ = &(selected_.back());
      }
      else if(gp->pdgId() == -6) {
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
      else if(Abs(gp->pdgId()) < 6 && Abs(gps[gp->momIdx()[0]].pdgId()) == 24) {
        selected_.push_back(*gp);
        wpartons_.push_back(&(selected_.back()));
      }
    }

    if(gp->momIdx().size() != 0 && Abs(gps[gp->momIdx()[0]].pdgId()) == 24) {	
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

bool  TTGenParticleSelector::select(URStreamer& event)
{
  //RESETS NEEDED VECTORS
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

  //SELECT BASED ON SELECTION MODE
  //TODO
  switch(mode_) {
  case NORMAL: select_normal(event); break;
  case PSEUDOTOP: select_pstop(event); break;
  case HERWIGPP: select_herwig(event); break;
    //case FULLDEP: select_with_deps(event); break;
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

/*
int TTGenParticleSelector::Collapse(int root, std::vector<const Genparticle*> &particles)
{
	bool found = true;
	while(found){
		found = false;
		for(auto particle : particles){
			if(particle->momIdx().size() == 0) continue;
			if(particle->momIdx()[0] == root){
				root = particle->idx();
				found = true;
				break;
			}
		}
	}
	return root;
}


void TTGenParticleSelector::select_with_deps(URStreamer& event)
{
	const std::vector<Genparticle>& gps = event.genParticles();
	TTbarHypothesis ret;

	//find top and tbar (multiple due to radiation
	std::vector<const Genparticle*> top_ids;
	std::vector<const Genparticle*> tbar_ids;
	std::vector<int> root_t, root_tbar;
	for(auto gp = gps.begin(); gp != gps.end(); ++gp){
		if(gp->pdgId() == ura::PDGID::t){
			top_ids.push_back(&(*gp));
			if(gp->momIdx().size() == 0 || gps[gp->momIdx()[0]].pdgId() != ura::PDGID::t){
				root_t.push_back(gp->idx());
			}
		} 
		else if(gp->pdgId() == ura::PDGID::tbar) {
			tbar_ids.push_back(&(*gp));
			if(gp->momIdx().size() == 0 || gps[gp->momIdx()[0]].pdgId() != ura::PDGID::tbar){
				root_tbar.push_back(gp->idx());
			}
		}
	}
	std::vector<int> collapsed_t, collapsed_tbar;
	for(auto idx : root_t){
		collapsed_t.push_back(
			Collapse(idx, top_ids)
			);
	}
	for(auto idx : root_tbar){
		collapsed_tbar.push_back(
			Collapse(idx, tbar_ids)
			);
	}
	//the std::vectors should be already sorted
	if(!is_sorted(collapsed_t.begin(), collapsed_t.end())) sort(collapsed_t.begin(), collapsed_t.end());
	if(!is_sorted(collapsed_tbar.begin(), collapsed_tbar.end())) sort(collapsed_tbar.begin(), collapsed_tbar.end());
		
	//look for top decay products
	int nbs = 0;
	int nbbars = 0;
	int nwp = 0, nwm = 0;
	int wpIdx = -1, wmIdx = -1;
	std::vector<const Genparticle*> wplus, wminus;
	for(auto gp = gps.begin(); gp != gps.end(); ++gp){
		if(gp->momIdx().size() == 0) continue;
		if(gp->pdgId() == ura::PDGID::b && binary_search(collapsed_t.begin(), collapsed_t.end(), gp->momIdx()[0])) {
			nbs++; 
			ret.b = &(*gp);
		}
		else if(gp->pdgId() == ura::PDGID::bbar && binary_search(collapsed_tbar.begin(), collapsed_tbar.end(), gp->momIdx()[0])) {
			nbbars++;
			ret.bbar = &(*gp);
		}
		else if(gp->pdgId() == ura::PDGID::Wminus){
			wminus.push_back(&(*gp));
			if(binary_search(collapsed_tbar.begin(), collapsed_tbar.end(), gp->momIdx()[0])) {nwm++; wmIdx = gp->idx();}
		}
		else if(gp->pdgId() == ura::PDGID::Wplus){
			wplus.push_back(&(*gp));
			if(binary_search(collapsed_t.begin(), collapsed_t.end(), gp->momIdx()[0])) {nwp++; wpIdx = gp->idx();}
		}
	}
	if(nbs != 1 && nbbars != 1)	Logger::log().error() << event.run<<":"<< event.lumi << ":" << event.evt << 
																" Found " << nbs << " b's and " << nbbars << " bbar's" << endl;
	if(nwp != 1 && nwm != 1) Logger::log().error() << event.run<<":"<< event.lumi << ":" << event.evt << 
														 " Found " << nwp << " W+'s and " << nwm << " W-'s" << endl;

	//collapse W+
	wpIdx = Collapse(wpIdx, wplus);
	wmIdx = Collapse(wmIdx, wminus);

	int wp_nprods = 0, wm_nprods = 0;
	const Genparticle *wplus_prods[2], *wminus_prods[2];		
	for(auto gp = gps.begin(); gp != gps.end(); ++gp){
		if(gp->momIdx().size() == 0) continue;
		if(gp->momIdx()[0] == wpIdx) {
			if(wp_nprods < 2) wplus_prods[wp_nprods] = &(*gp);
			wp_nprods++; 
		}
		else if(gp->momIdx()[0] == wmIdx) {
			if(wm_nprods < 2) wminus_prods[wm_nprods] = &(*gp);
			wm_nprods++; 
		}
	}
	if(wp_nprods != 2 && wm_nprods != 2){ 
		Logger::log().error() << event.run<<":"<< event.lumi << ":" << event.evt << 
			" Found " << wp_nprods <<
			" W+ products and " << wm_nprods << " W- products" << endl;
		ret.decay = INVALID;
	} else {
		int lep_decays = 0;
		//w products are quarks
		bool wp1_isquark = fabs(wplus_prods[0]->pdgId()) <= 6;
		bool wp2_isquark = fabs(wplus_prods[1]->pdgId()) <= 6;
		//w products are leptons or quarks
		bool wp1_isqlep = fabs(wplus_prods[0]->pdgId()) <= 16;
		bool wp2_isqlep = fabs(wplus_prods[1]->pdgId()) <= 16;
		if(wp1_isquark && wp2_isquark) {
			ret.wplus.first  = wplus_prods[0];
			ret.wplus.second = wplus_prods[1];
			ret.wplus.isLeptonic = false;
		} else if((wp1_isqlep && !wp1_isquark) && (wp2_isqlep && !wp2_isquark)) {
			//assign the charged lepton as first
			int id0 = wplus_prods[0]->pdgId();
			ret.wplus.first  = (id0 % 2 == 0) ? wplus_prods[1] : wplus_prods[0];
			ret.wplus.second = (id0 % 2 == 0) ? wplus_prods[0] : wplus_prods[1];
			lep_decays++;
			ret.wplus.isLeptonic = true;
		} else {
			Logger::log().error() << event.run<<":"<< event.lumi << ":" << event.evt << 
				" W+ decays to lepton and quark! (" <<
				wplus_prods[0]->pdgId() << ", " << wplus_prods[1]->pdgId() << endl;
			// ") 1:" << wp1_isquark << "," << wp1_isqlep << " 2:" <<
			// wp1_isquark << "," << wp2_isqlep << endl;
			ret.decay = INVALID;
		}

		//FIXME: this should go in a separate function
		//w products are quarks
		bool wm1_isquark = fabs(wminus_prods[0]->pdgId()) <= 6;
		bool wm2_isquark = fabs(wminus_prods[1]->pdgId()) <= 6;
		//w products are leptons or quarks
		bool wm1_isqlep = fabs(wminus_prods[0]->pdgId()) <= 16;
		bool wm2_isqlep = fabs(wminus_prods[1]->pdgId()) <= 16;
		if(wm1_isquark && wm2_isquark) {
			ret.wminus.first  = wminus_prods[0];
			ret.wminus.second = wminus_prods[1];
			ret.wminus.isLeptonic = false;
		} else if((wm1_isqlep && !wm1_isquark) && (wm2_isqlep && !wm2_isquark)) {
			//assign the charged lepton as first
			int id0 = wminus_prods[0]->pdgId();
			ret.wminus.first  = (id0 % 2 == 0) ? wminus_prods[1] : wminus_prods[0];
			ret.wminus.second = (id0 % 2 == 0) ? wminus_prods[0] : wminus_prods[1];
			lep_decays++;
			ret.wminus.isLeptonic = true;
		} else {
			Logger::log().error() << event.run<<":"<< event.lumi << ":" << event.evt << 
				" W- decays to lepton and quark! (" <<
				wminus_prods[0]->pdgId() << ", " << wminus_prods[1]->pdgId() <<
				")" << endl;
			ret.decay = INVALID;
		}
		if(ret.decay != INVALID){
			switch(lep_decays){
			case 0: ret.decay = FULLHAD; break;
			case 1: ret.decay = SEMILEP; break;
			case 2: ret.decay = FULLLEP; break;
			}
		}
	}
	return ret;

}//*/

