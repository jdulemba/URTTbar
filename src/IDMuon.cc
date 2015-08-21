#include <IDMuon.h>
#include <TMath.h>
#include <iostream>
#include "Logger.h"

using namespace TMath;

IDMuon::IDMuon(const Muon mu, double rho): 
	Muon(mu), 
	MCMatchable(),
	rho_(rho)
{
}

double IDMuon::PFIsoDb()
{
	return (chargedIso() + TMath::Max(neutralIso() + photonIso() - 0.5*puIso(), 0.));
}
	
double IDMuon::CorPFIsolation2015()
{
	double eta = Abs(Eta());
	double effarea = 0.;
	if(eta < 0.8){ effarea = 0.0913;}
	else if(eta < 1.3){ effarea = 0.0765;}
	else if(eta < 2.0){ effarea = 0.0546;}
	else if(eta < 2.2){ effarea = 0.0728;}
	else if(eta < 2.5){ effarea = 0.1177;}

	if(rho_ < 0.){Logger::log().error() << "Store the value of rho in the electrons to use this isolation" << endl;}
	effarea *= Max(rho_, 0.);
	return(chargedIso() + Max(neutralIso() + photonIso() - effarea, 0.))/Pt();
}

bool IDMuon::ID(IDS idtyp)
{
	if(idtyp == TIGHT_12 || idtyp == TIGHT_12Db)
	{
		if(TMath::Abs(Eta()) > 2.4) return(false);
		if(!isPF()) return(false);
		if(!isGlobal()) return(false);
		if(pixelHits() <= 0) return(false);
		if(trackerLayers() <= 5) return(false);
		if(validHits() <= 0) return(false);
		if(numMatchedStations() <= 1) return(false);
		if(TMath::Abs(dB()) > 0.2) return(false);
		if(TMath::Abs(dz()) > 0.5) return(false);
		if(chi2()/ndof() > 10.) return(false);
		if(USEISO && idtyp == TIGHT_12Db && PFIsoDb()/Pt() > 0.12) return(false);
		if(USEISO && idtyp == TIGHT_12 && (trackiso())/Pt() > 0.05) return(false);
		//if(idtyp == TIGHT_12 && CorPFIsolation2015()/Pt() > 0.05) return(false);
		return(true);
	}
	else if(idtyp == LOOSE_12 || idtyp == LOOSE_12Db)
	{
		if(TMath::Abs(Eta()) > 2.4) return(false);
		if(!isPF()) return(false);
		if(!isGlobal() && !isTracker()) return(false);
		if(USEISO && idtyp == LOOSE_12Db && PFIsoDb()/Pt() > 0.2) return(false);
		if(USEISO && idtyp == LOOSE_12 && (trackiso())/Pt() > 0.1) return(false);
		//if(idtyp == LOOSE_12 && CorPFIsolation2015()/Pt() > 0.1) return(false);
		return(true);
	}
	return(false);
}

bool IDMuon::USEISO = true;
