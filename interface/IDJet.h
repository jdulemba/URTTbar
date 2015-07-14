#ifndef IDJET_H
#define IDJET_H
#include "MCMatchable.h"
#include "URStreamer.h"
#include "IDMuon.h"
#include "IDElectron.h"
#include <TMath.h>

class IDJet : public Jet, public MCMatchable
{
private:
	double rndm_;
public:
	enum BTag {CSVLOOSE, CSVMEDIUM, CSVTIGHT};
	IDJet(const Jet el, double rndm):
		Jet(el),
    MCMatchable(),
		rndm_(rndm)
		{
		}

	IDJet(const Jet el):
		Jet(el),
    MCMatchable(),
		rndm_(-1)
		{
		}

	int flavor() const {return (match()) ? match()->pdgId() : partonFlavour();}

	double rndm() const {return rndm_;}

	bool BTagId(BTag wp) const
	{
		double threshold = -1.;
		if(wp == BTag::CSVLOOSE) threshold = 0.605;
		else if(wp == BTag::CSVMEDIUM) threshold = 0.890;
		else if(wp == BTag::CSVTIGHT) threshold = 0.970;
		
		return csvIncl() > threshold;
	}

	bool ID()
	{
		//to be filled in new tree version
		if(numberOfDaughters() <= 1) {return false;}
		if(neutralHadronEnergyFraction() + HFHadronEnergyFraction() >= 0.99){return false;}
		if(neutralEmEnergyFraction() >= 0.99){return false;}
		if(TMath::Abs(Eta()) < 2.4)
		{
			if(chargedEmEnergyFraction() >= 0.99){return false;}
			if(chargedHadronEnergyFraction() <= 0.){return false;}
			if(chargedMultiplicity() <= 0.){return false;}
		}
		return(true);
	}
	bool Clean(const vector<IDMuon*>& muons, const vector<IDElectron*>& electrons, double distpar = 0.4)
	{
		for(const IDMuon* mu : muons)
		{
			if(DeltaR(*mu) < distpar)
			{
				return false;
			}
		}
		for(const IDElectron* el : electrons)
		{
			if(DeltaR(*el) < distpar)
			{
				return false;
			}
		}
		return true;
	}

};
#endif
