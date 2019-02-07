#include "Analyses/URTTbar/interface/IDElectron.h"
#include <TMath.h>
#include <iostream>
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include <cmath>
using namespace TMath;

const std::map<std::string, IDElectron::IDS> IDElectron::id_names = {
    {"FAIL"   , IDElectron::IDS::FAIL   },
    {"TIGHT_15"   , IDElectron::IDS::TIGHT_15   },
    {"MEDIUM_15"  , IDElectron::IDS::MEDIUM_15  },
    {"LOOSE_15"   , IDElectron::IDS::LOOSE_15   },
    {"VETO_15"    , IDElectron::IDS::VETO_15    },
    {"TIGHT_15_NoECAL_Gap", IDElectron::IDS::TIGHT_15_NoECAL_Gap},
    {"NOTVETO_15", IDElectron::IDS::NOTVETO_15},
    {"FAKES", IDElectron::IDS::FAKES},
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
    //EA Values available here: https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
    double abseta = fabs(eta);
    if(abseta >= 0.0000 && abseta <= 1.0000) effarea = 0.1703;
    if(abseta >  1.0000 && abseta <= 1.4790) effarea = 0.1715;
    if(abseta >  1.4790 && abseta <= 2.0000) effarea = 0.1213;
    if(abseta >  2.0000 && abseta <= 2.2000) effarea = 0.1230;
    if(abseta >  2.2000 && abseta <= 2.3000) effarea = 0.1635;
    if(abseta >  2.3000 && abseta <= 2.4000) effarea = 0.1937;
    if(abseta >  2.4000 && abseta <= 5.0000) effarea = 0.2393;

    if(rho_ < 0.){Logger::log().error() << "Store the value of rho in the electrons to use this isolation: " << rho_ << endl;}
    effarea *= Max(rho_, 0.);
    //return((PFR3().Charged() + Max(PFR3().Neutral() + PFR3().Photon() - Max(GLAN->AK5PFRho(), 0.f)*effarea, 0.))/Pt());
    //return(chargedIso() + Max(neutralIso() + photonIso() - effarea, 0.))/Pt();
    return(pfHadronIso() + Max(pfNeutralIso() + pfPhotonIso() - effarea, 0.));
}

bool IDElectron::FakeID() const {// {return IPCuts() && (eidCutNoIsoTight() > 0.5) ;TIGHT_15_NoECAL_Gap
    bool ecalgap = (fabs(etaSC()) <= 1.4442 || fabs(etaSC()) >= 1.5660);
    bool iso = (etaSC() < 1.479) ? PFIsolationRho2015() >= 0.0588 : PFIsolationRho2015() >= 0.0571;
    return IPCuts() && eidCutNoIsoTight() && ecalgap && iso;
}

bool IDElectron::ID(IDS idtyp)
{
    double sceta = Abs(TVector3(x(), y(), z()).Eta());
    if(sceta > 2.5) return(false);
    switch(idtyp) {
        case FAIL: return false;
        case TIGHT_15: return TightID25ns();
        case MEDIUM_15: return MediumID25ns();
        case LOOSE_15: return LooseID25ns();
        case VETO_15: return VetoID25ns();
        case TIGHT_15_NoECAL_Gap: return (TightID25ns() && (fabs(etaSC()) <= 1.4442 || fabs(etaSC()) >= 1.5660)); //removes EB-EE gap  
        case NOTVETO_15: return !VetoID25ns();
        case FAKES: return FakeID();

    }
    return false;
}

bool IDElectron::USEISO = true;

