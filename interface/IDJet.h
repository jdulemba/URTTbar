#ifndef IDJET_H
#define IDJET_H
#include "MCMatchable.h"
#include "URStreamer.h"
#include "IDMuon.h"
#include "IDElectron.h"
#include <TMath.h>
#include <unordered_map>
#include <string>
#include "Logger.h"
#include "BTagCalibrationStandalone.h"

class IDJet : public Jet, public MCMatchable
{
private:
	double rndm_;
  TLorentzVector uncorr_;
public:
	enum BTag {NONE, CSVLOOSE, CSVMEDIUM, CSVTIGHT, CTAGLOOSE, CTAGMEDIUM, CTAGTIGHT};
	enum IDType {NOTSET, CSV, CTAG};
  
  static const std::unordered_map<std::string, IDJet::BTag> tag_names;
  static IDJet::BTag tag(std::string label);

	IDJet(const Jet el, double rndm):
		Jet(el),
    MCMatchable(),
		rndm_(rndm),
    uncorr_(el)
		{
		}

	IDJet(const Jet el):
		Jet(el),
    MCMatchable(),
		rndm_(-1),
    uncorr_(el)
		{
		}

	int flavor() const {return (match()) ? match()->pdgId() : partonFlavour();}
  void update_energy(float val) {SetPtEtaPhiE(val*TMath::Sin(Theta()), Eta(), Phi(), val);}
  void resetp4() {SetPtEtaPhiE(uncorr_.Pt(), uncorr_.Eta(), uncorr_.Phi(), uncorr_.E());}
  TLorentzVector original_p4() {return uncorr_;}

	double rndm() const {return rndm_;}

  static std::string tag2string(BTag id);
  static IDType id_type(BTag id);
  static std::string id_string(BTag id);
  static BTagEntry::OperatingPoint tag_tightness(BTag id);

	bool BTagId(BTag wp) const;
	bool CTagId(BTag wp) const;
  bool TagId(BTag wp) const;

	bool ID()	{
		//to be filled in new tree version
		if(numberOfDaughters() <= 1) {return false;}
		if(neutralHadronEnergyFraction() + HFHadronEnergyFraction() >= 0.99){return false;}
		if(neutralEmEnergyFraction() >= 0.99){return false;}
		if(TMath::Abs(Eta()) < 2.4)	{
			if(chargedEmEnergyFraction() >= 0.99){return false;}
			if(chargedHadronEnergyFraction() <= 0.){return false;}
			if(chargedMultiplicity() <= 0.){return false;}
		}
		return(true);
	}

	bool Clean(const vector<IDMuon*>& muons, const vector<IDElectron*>& electrons, double distpar = 0.4) {
		for(const IDMuon* mu : muons) {
			if(DeltaR(*mu) < distpar) {
				return false;
			}
		}
		for(const IDElectron* el : electrons)	{
			if(DeltaR(*el) < distpar)	{
				return false;
			}
		}
		return true;
	}

};

#endif
