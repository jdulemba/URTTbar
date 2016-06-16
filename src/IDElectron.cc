#include "Analyses/URTTbar/interface/IDElectron.h"
#include <TMath.h>
#include <iostream>
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include <cmath>
using namespace TMath;

const std::map<std::string, IDElectron::IDS> IDElectron::id_names = {
  {"TIGHT_15"   , IDElectron::IDS::TIGHT_15   },
  {"MEDIUM_15"  , IDElectron::IDS::MEDIUM_15  },
  {"LOOSE_15"   , IDElectron::IDS::LOOSE_15   }
};

IDElectron::IDS IDElectron::id(const std::string label) {
  try {
    return IDElectron::id_names.at(label);
  }
  catch (std::out_of_range e){
    Logger::log().error() << "Requested IDElectron ID label: " << label << " does not exist!" << std::endl;
    throw 49;
  }
}

double IDElectron::PFIsolationRho2015() const
{
  // Isolation is calculated following an example in [1], as recommended in [2]
	//[1] https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_2/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleEffAreaPFIsoCut.cc#L75-L90
	//[2] https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1664/1.html
	double eta = etaSC(); //Abs(TVector3(x(), y(), z()).Eta());
	double effarea = 0.;
  // The following values refer to EA for cone 0.3 and fixedGridRhoFastjetAll. 
  // They are valid for electrons only, different EA are available for muons.
	//EA Values available here: https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt
  if(abs(eta)>=0.0 && abs(eta)<=1.0) effarea = 0.1752;
  if(abs(eta)>1.0 && abs(eta)<=1.479) effarea = 0.1862;
  if(abs(eta)>1.479 && abs(eta)<=2.0) effarea = 0.1411;
  if(abs(eta)>2.0 && abs(eta)<=2.2) effarea = 0.1534;
  if(abs(eta)>2.2 && abs(eta)<=2.3) effarea = 0.1903;
  if(abs(eta)>2.3 && abs(eta)<=2.4) effarea = 0.2243;
  if(abs(eta)>2.4 && abs(eta)<=5 ) effarea = 0.2687;

	if(rho_ < 0.){Logger::log().error() << "Store the value of rho in the electrons to use this isolation" << endl;}
	effarea *= Max(rho_, 0.);
	//return((PFR3().Charged() + Max(PFR3().Neutral() + PFR3().Photon() - Max(GLAN->AK5PFRho(), 0.f)*effarea, 0.))/Pt());
	//return(chargedIso() + Max(neutralIso() + photonIso() - effarea, 0.))/Pt();
	return(pfHadronIso() + Max(pfNeutralIso() + pfPhotonIso() - effarea, 0.));
}

//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
bool IDElectron::LooseID25ns() const {
  return (eidCutLoose() > 0.5);
  //if(full5x5SigmaIEtaIEta() > 0.01){return(false);}
  // if(full5x5_sigmaIEtaIEta()  >= ((isEB()) ? 0.0103: 0.0301)) return false; //to be updated to full5x5SigmaIEtaIEta
  // if(Abs(DEtaSCTrk()) >= ((isEB()) ? 0.0105: 0.00814)) return false;
  // if(Abs(DPhiSCTrk()) >= ((isEB()) ? 0.115:  0.182)) return false;
  // if(hadronicOverEM() >= ((isEB()) ? 0.104:  0.0897)) return false;
  // if(PFIsolationRho2015() >= ((isEB()) ? 0.0893: 0.121)) return false; //to be updated
  // //relIsoWithEA
  // float ooEmooP = (ecalEnergy() == 0 || !std::isfinite(ecalEnergy())) ? 999 : Abs(1.0/ecalEnergy() - ESCOverETrack()/ecalEnergy() );
  // if(ooEmooP    >= ((isEB()) ? 0.102 : 0.126)) return false;
  // if(Abs(dz())  >= ((isEB()) ? 0.0261: 0.118)) return false; //to be checked, needs vtx?
  // if(Abs(dB()) >= ((isEB()) ? 0.41  : 0.822)) return false; //to be checked, needs vtx?
  // if(nMissingInnerHits() > ((isEB()) ? 2: 1)) return false; //to be updated to nMissingTrackerHits
  // if(!passConversionVeto()) return false;
  // return true;
}

bool IDElectron::MediumID25ns() const {
  return (eidCutMedium() > 0.5);
  //if(full5x5SigmaIEtaIEta() > 0.01){return(false);}
  // if(sigmaIEtaIEta()  >= ((isEB()) ? 0.0101: 0.0283 )) return false; //to be updated to full5x5SigmaIEtaIEta
  // if(Abs(DEtaSCTrk()) >= ((isEB()) ? 0.0103: 0.00733)) return false;
  // if(Abs(DPhiSCTrk()) >= ((isEB()) ? 0.0336: 0.114  )) return false;
  // if(hadronicOverEM() >= ((isEB()) ? 0.0876: 0.0678 )) return false;
  // if(PFIsolationRho2015() >= ((isEB()) ? 0.0766: 0.0678 )) return false; //to be updated
  // //relIsoWithEA
  // float ooEmooP = (ecalEnergy() == 0 || !std::isfinite(ecalEnergy())) ? 999 : Abs(1.0/ecalEnergy() - ESCOverETrack()/ecalEnergy() );
  // if(ooEmooP    >= ((isEB()) ? 0.0174: 0.0898)) return false;
  // if(Abs(dz())  >= ((isEB()) ? 0.0118: 0.0739)) return false; //to be checked, needs vtx?
  // if(Abs(dB()) >= ((isEB()) ? 0.373 : 0.602 )) return false; //to be checked, needs vtx?
  // if(nMissingInnerHits() > ((isEB()) ? 2: 1)) return false; //to be updated to nMissingTrackerHits
  // if(!passConversionVeto()) return false;
  // return true;
}

bool IDElectron::TightID25ns() const {
  return (eidCutTight() > 0.5);
  //if(full5x5SigmaIEtaIEta() > 0.01){return(false);}
  // if(sigmaIEtaIEta()  >= ((isEB()) ? 0.0101 : 0.0279 )) return false; //to be updated to full5x5SigmaIEtaIEta
  // if(Abs(DEtaSCTrk()) >= ((isEB()) ? 0.00926: 0.00724)) return false;
  // if(Abs(DPhiSCTrk()) >= ((isEB()) ? 0.0336 : 0.0918 )) return false;
  // if(hadronicOverEM() >= ((isEB()) ? 0.0597 : 0.0615 )) return false;
  // if(PFIsolationRho2015() >= ((isEB()) ? 0.0354 : 0.0646 )) return false; //to be updated
  // //relIsoWithEA
  // float ooEmooP = (ecalEnergy() == 0 || !std::isfinite(ecalEnergy())) ? 999 : Abs(1.0/ecalEnergy() - ESCOverETrack()/ecalEnergy() );
  // if(ooEmooP    >= ((isEB()) ? 0.012 : 0.00999)) return false;
  // if(Abs(dz())  >= ((isEB()) ? 0.0111: 0.0351 )) return false; //to be checked, needs vtx?
  // if(Abs(dB()) >= ((isEB()) ? 0.0466: 0.417  )) return false; //to be checked, needs vtx?
  // if(nMissingInnerHits() > ((isEB()) ? 2: 1)) return false; //to be updated to nMissingTrackerHits
  // if(!passConversionVeto()) return false;
  // return true;
}

bool IDElectron::ID(IDS idtyp)
{
  double sceta = Abs(TVector3(x(), y(), z()).Eta());
	if(sceta > 2.5) return(false);
  switch(idtyp) {
  case TIGHT_15: return TightID25ns();
  case MEDIUM_15: return MediumID25ns();
  case LOOSE_15: return LooseID25ns();
  }
  return false;
}

bool IDElectron::USEISO = true;

